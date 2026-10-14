#include "gcrank/gc_rank.h"

#include <algorithm>
#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <mutex>
#include <print>
#include <span>
#include <string_view>
#include <thread>
#include <utility>

#include "gcrank/crc32.h"
#include "gcrank/error.h"
#include "gcrank/fasta_reader.h"
#include "gcrank/format.h"
#include "gcrank/input_stream.h"
#include "gcrank/sha256.h"
#include "gcrank/thread_pool.h"

namespace gcrank {

namespace fs = std::filesystem;

unsigned default_threads() {
  const unsigned hardware = std::thread::hardware_concurrency();
  const unsigned available = hardware == 0 ? 4u : hardware;
  return std::min(8u, std::max(1u, available));
}

namespace {

// ---- parallel build pipeline --------------------------------------------
//
//   producer (1 thread)  parses FASTA -> bounded ring of slots (chunk_bits bases)
//   workers (N threads)  pack A/C/G/T/N bases into GC and valid words
//   assembler (caller)   appends blocks in order, records contig boundaries
//
// Slots are recycled through the ring; the assembler consumes in submission
// order, so the output does not depend on the thread count.

struct Slot {
  enum class Kind { ContigBegin, ContigEnd, Block };

  Kind kind = Kind::Block;
  std::string name;            // ContigBegin
  uint64_t nbases = 0;         // Block
  std::vector<char> bases;     // Block, producer filled
  std::vector<uint64_t> gc_words;     // Block, worker filled
  std::vector<uint64_t> valid_words;  // Block, worker filled

  bool free = true;   // producer may fill it
  bool ready = false; // assembler may consume it
};

// Packs bases into MSB-first 64-bit words (bit 0 of the block = MSB of word 0).
// Padding bits after `n` bases in the last word are zero.
void pack_block(std::span<const char> bases, std::vector<uint64_t>& gc_words,
                std::vector<uint64_t>& valid_words) {
  const std::size_t n = bases.size();
  const std::size_t nwords = (n + 63) / 64;
  gc_words.assign(nwords, 0);
  valid_words.assign(nwords, 0);
  std::size_t word = 0;
  uint64_t gc_acc = 0;
  uint64_t valid_acc = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const BaseClass cls = base_class(bases[i]);
    gc_acc = (gc_acc << 1) | (base_is_gc(cls) ? uint64_t(1) : uint64_t(0));
    valid_acc = (valid_acc << 1) | (base_is_valid(cls) ? uint64_t(1) : uint64_t(0));
    if ((i & 63u) == 63u) {
      gc_words[word] = gc_acc;
      valid_words[word] = valid_acc;
      ++word;
      gc_acc = 0;
      valid_acc = 0;
    }
  }
  const unsigned rem = static_cast<unsigned>(n & 63u);
  if (rem != 0) {
    gc_words[word] = gc_acc << (64 - rem);
    valid_words[word] = valid_acc << (64 - rem);
  }
}

// Thrown inside the producer thread when the assembler is unwinding.
struct PipelineAbort {};

class BuildPipeline {
public:
  BuildPipeline(std::string path, uint64_t chunk_bases, unsigned workers, bool verbose)
      : path_(std::move(path)),
        chunk_(chunk_bases),
        workers_(workers),
        verbose_(verbose),
        capacity_(static_cast<uint64_t>(workers) + 2),
        slots_(static_cast<std::size_t>(capacity_)),
        pool_(workers) {}

  void run(BitVector& genome, BitVector& valid, std::vector<ContigInfo>& contigs) {
    genome_ = &genome;
    valid_ = &valid;
    contigs_ = &contigs;
    producer_ = std::jthread([this](std::stop_token stop) { produce(stop); });
    for (unsigned i = 0; i < workers_; ++i) pool_.submit([this] { work_loop(); });
    try {
      assemble();
    } catch (...) {
      shutdown();
      throw;
    }
    shutdown();
  }

private:
  // ---- producer ----------------------------------------------------------

  void produce(std::stop_token stop) {
    try {
      std::unique_ptr<FastaReader> reader = FastaReader::from_path(path_, chunk_);
      uint64_t index = 0;
      uint64_t bases = 0;
      auto next_progress = std::chrono::steady_clock::now() + std::chrono::seconds(4);
      for (;;) {
        const FastaReader::Event event = reader->next();
        if (event == FastaReader::Event::End) break;

        Slot& slot = slots_[static_cast<std::size_t>(index % capacity_)];
        {
          std::unique_lock<std::mutex> lock(mutex_);
          cv_free_.wait(lock, [&] {
            return slot.free || failed_.load() || stop.stop_requested();
          });
          if (failed_.load() || stop.stop_requested()) throw PipelineAbort{};
          slot.free = false;
        }

        bool is_block = false;
        switch (event) {
          case FastaReader::Event::ContigBegin:
            slot.kind = Slot::Kind::ContigBegin;
            slot.name = reader->contig_name();
            break;
          case FastaReader::Event::Block:
            slot.kind = Slot::Kind::Block;
            slot.bases.assign(reader->block().begin(), reader->block().end());
            slot.nbases = static_cast<uint64_t>(slot.bases.size());
            bases += slot.nbases;
            is_block = true;
            break;
          case FastaReader::Event::ContigEnd:
            slot.kind = Slot::Kind::ContigEnd;
            slot.name.clear();
            break;
          default:
            break;
        }
        publish(index, is_block);
        ++index;

        if (verbose_ && is_block) {
          const auto now = std::chrono::steady_clock::now();
          if (now >= next_progress) {
            next_progress = now + std::chrono::seconds(4);
            std::lock_guard<std::mutex> log(verbose_mutex_);
            std::print(stderr, "[gcrank] parsed {} bases\r", bases);
            std::fflush(stderr);
          }
        }
      }
      if (verbose_) {
        std::lock_guard<std::mutex> log(verbose_mutex_);
        std::print(stderr, "[gcrank] parsed {} bases\n", bases);
        std::fflush(stderr);
      }
    } catch (const PipelineAbort&) {
      // The assembler is failing; unwind quietly.
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!producer_error_ && !failed_.load()) producer_error_ = std::current_exception();
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      producer_done_ = true;
    }
    cv_ready_.notify_all();
  }

  void publish(uint64_t index, bool is_block) {
    if (is_block) {
      {
        std::lock_guard<std::mutex> lock(mutex_);
        produce_index_ = index + 1;
      }
      {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        work_.push_back(index);
      }
      cv_work_.notify_one();
    } else {
      {
        std::lock_guard<std::mutex> lock(mutex_);
        produce_index_ = index + 1;
        slots_[static_cast<std::size_t>(index % capacity_)].ready = true;
      }
      cv_ready_.notify_one();
    }
  }

  // ---- workers -----------------------------------------------------------

  void work_loop() {
    for (;;) {
      uint64_t index = 0;
      {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        cv_work_.wait(lock, [&] { return !work_.empty() || workers_stop_ || failed_.load(); });
        if (work_.empty()) return;
        index = work_.front();
        work_.pop_front();
      }
      Slot& slot = slots_[static_cast<std::size_t>(index % capacity_)];
      try {
        pack_block(slot.bases, slot.gc_words, slot.valid_words);
        // Keep the capacity for reuse; this is the buffer the peak-memory
        // budget is about.
        slot.bases.clear();
      } catch (...) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (!worker_error_) worker_error_ = std::current_exception();
        failed_ = true;
        cv_ready_.notify_all();
        return;
      }
      {
        std::lock_guard<std::mutex> lock(mutex_);
        slot.ready = true;
      }
      cv_ready_.notify_one();
    }
  }

  // ---- assembler (calling thread) ----------------------------------------

  void assemble() {
    ContigInfo current;
    for (;;) {
      std::exception_ptr error;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_ready_.wait(lock, [&] {
          if (failed_.load() || producer_error_) return true;
          if (consume_index_ < produce_index_ &&
              slots_[static_cast<std::size_t>(consume_index_ % capacity_)].ready)
            return true;
          return producer_done_ && consume_index_ == produce_index_;
        });
        if (worker_error_) {
          error = worker_error_;
        } else if (producer_error_) {
          error = producer_error_;
        } else if (failed_.load()) {
          error = std::make_exception_ptr(Error("gcrank: internal pipeline failure"));
        } else if (producer_done_ && consume_index_ == produce_index_) {
          break;  // all items consumed
        }
      }
      if (error) std::rethrow_exception(error);

      Slot& slot = slots_[static_cast<std::size_t>(consume_index_ % capacity_)];
      switch (slot.kind) {
        case Slot::Kind::ContigBegin:
          current.name = slot.name;
          current.length = 0;
          current.offset = genome_->size();
          break;
        case Slot::Kind::Block:
          genome_->append_bits(slot.gc_words, slot.nbases);
          valid_->append_bits(slot.valid_words, slot.nbases);
          current.length += slot.nbases;
          break;
        case Slot::Kind::ContigEnd:
          contigs_->push_back(current);
          break;
      }
      {
        std::lock_guard<std::mutex> lock(mutex_);
        slot.ready = false;
        slot.free = true;
        ++consume_index_;
      }
      cv_free_.notify_one();
    }
  }

  void shutdown() {
    producer_.request_stop();
    {
      std::lock_guard<std::mutex> lock(queue_mutex_);
      workers_stop_ = true;
    }
    cv_work_.notify_all();
    cv_free_.notify_all();
    cv_ready_.notify_all();
    if (producer_.joinable()) producer_.join();
    try {
      pool_.wait();
      pool_.stop();
    } catch (...) {
      // Worker exceptions are already captured in worker_error_.
    }
  }

  std::string path_;
  uint64_t chunk_;
  unsigned workers_;
  bool verbose_;
  uint64_t capacity_;
  std::vector<Slot> slots_;

  std::mutex mutex_;
  std::condition_variable cv_free_;
  std::condition_variable cv_ready_;
  uint64_t produce_index_ = 0;
  uint64_t consume_index_ = 0;
  bool producer_done_ = false;
  std::exception_ptr producer_error_;
  std::exception_ptr worker_error_;
  std::atomic<bool> failed_{false};

  std::mutex queue_mutex_;
  std::condition_variable cv_work_;
  std::deque<uint64_t> work_;
  bool workers_stop_ = false;

  std::mutex verbose_mutex_;
  std::jthread producer_;
  ThreadPool pool_;

  BitVector* genome_ = nullptr;
  BitVector* valid_ = nullptr;
  std::vector<ContigInfo>* contigs_ = nullptr;
};

// ---- index files ---------------------------------------------------------

FileInfo write_contigs(const fs::path& path, const std::vector<ContigInfo>& contigs) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path.string() + "' for writing");
  CRC32 crc;
  uint64_t bytes = 0;
  std::string line;
  for (const ContigInfo& contig : contigs) {
    if (contig.name.find('\n') != std::string::npos)
      throw Error("contig name contains a newline; refusing to write '" + path.string() + "'");
    line.clear();
    std::format_to(std::back_inserter(line), "{}\t{}\t{}\n", contig.name, contig.length,
                   contig.offset);
    out.write(line.data(), static_cast<std::streamsize>(line.size()));
    if (!out) throw OutputError("write error on '" + path.string() + "'");
    crc.update(line.data(), line.size());
    bytes += line.size();
  }
  out.flush();
  if (!out) throw OutputError("write error on '" + path.string() + "'");
  return FileInfo{bytes, crc.value()};
}

uint64_t parse_u64_field(std::string_view text, const std::string& path, std::string_view what) {
  if (text.empty()) throw FormatError("'" + path + "': empty " + std::string(what));
  uint64_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec == std::errc::result_out_of_range)
    throw FormatError("'" + path + "': " + std::string(what) + " overflows: '" +
                      std::string(text) + "'");
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw FormatError("'" + path + "': " + std::string(what) + " is not a number: '" +
                      std::string(text) + "'");
  return value;
}

std::vector<ContigInfo> parse_contigs(const fs::path& path, uint64_t expected_contigs,
                                      uint64_t expected_bases) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path.string() + "' for reading");
  const std::string name = path.string();
  std::vector<ContigInfo> contigs;
  std::string line;
  uint64_t next_offset = 0;
  while (std::getline(in, line)) {
    std::string_view view(line);
    if (!view.empty() && view.back() == '\r') view.remove_suffix(1);
    if (view.empty()) continue;
    const std::size_t tab_offset = view.rfind('\t');
    if (tab_offset == std::string_view::npos)
      throw FormatError("'" + name + "': expected '<name>\\t<length>\\t<offset>'");
    const std::size_t tab_length = view.rfind('\t', tab_offset - 1);
    if (tab_length == std::string_view::npos)
      throw FormatError("'" + name + "': expected '<name>\\t<length>\\t<offset>'");
    ContigInfo contig;
    contig.name = view.substr(0, tab_length);
    contig.length = parse_u64_field(view.substr(tab_length + 1, tab_offset - tab_length - 1),
                                    name, "contig length");
    contig.offset = parse_u64_field(view.substr(tab_offset + 1), name, "contig offset");
    if (contig.offset != next_offset)
      throw FormatError("'" + name + "': contig '" + contig.name + "' starts at " +
                        std::to_string(contig.offset) + " but " +
                        std::to_string(next_offset) + " was expected");
    next_offset += contig.length;
    contigs.push_back(std::move(contig));
  }
  if (in.bad()) throw IoError("read error on '" + name + "'");
  if (contigs.size() != expected_contigs)
    throw FormatError("'" + name + "': holds " + std::to_string(contigs.size()) +
                      " contigs but the manifest declares " +
                      std::to_string(expected_contigs));
  if (next_offset != expected_bases)
    throw FormatError("'" + name + "': contigs cover " + std::to_string(next_offset) +
                      " bases but the manifest declares " + std::to_string(expected_bases));
  return contigs;
}

void write_sidecar(const fs::path& path, uint32_t crc) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path.string() + "' for writing");
  const std::string text = crc32_hex(crc) + "\n";
  out.write(text.data(), static_cast<std::streamsize>(text.size()));
  out.flush();
  if (!out) throw OutputError("write error on '" + path.string() + "'");
}

// Re-reads the data files and compares them against the manifest. This is what
// makes the temporary-directory -> final rename safe. Returns an empty string
// when everything matches.
std::string verify_index_files(const fs::path& dir, const Manifest& manifest) {
  for (const auto& entry : manifest.files) {
    const fs::path path = dir / entry.first;
    std::error_code ec;
    const uint64_t size = fs::file_size(path, ec);
    if (ec)
      return "cannot stat '" + path.string() + "': " + ec.message();
    if (size != entry.second.size)
      return "'" + path.string() + "' holds " + std::to_string(size) + " bytes but " +
             std::to_string(entry.second.size) + " were written";
    uint32_t crc = 0;
    try {
      crc = crc32_file(path.string());
    } catch (const IoError& error) {
      return error.what();
    }
    if (crc != entry.second.crc32)
      return "'" + path.string() + "' failed its checksum validation";
  }
  return std::string();
}

void load_bitvector(BitVector& bitvector, const fs::path& bv_path, const fs::path& rb_path,
                    uint64_t nbits, bool stored_little_endian, bool use_mmap) {
  if (use_mmap) {
    try {
      MappedFile bv_file = MappedFile::map(bv_path.string());
      MappedFile rb_file = MappedFile::map(rb_path.string());
      bitvector.map_bv(std::move(bv_file), nbits);
      bitvector.map_rank(std::move(rb_file), nbits, stored_little_endian);
      return;
    } catch (const IoError& error) {
      std::print(stderr, "gcrank: warning: {}; falling back to reading\n", error.what());
    }
  }
  bitvector.load_bv(bv_path.string(), nbits);
  bitvector.load_rank(rb_path.string(), nbits, stored_little_endian);
}

}  // namespace

// ---- GcRank ---------------------------------------------------------------

GcRank::GcRank() = default;
GcRank::~GcRank() = default;
GcRank::GcRank(GcRank&&) noexcept = default;
GcRank& GcRank::operator=(GcRank&&) noexcept = default;

void GcRank::build(const std::string& fasta_path, const BuildOptions& options) {
  if (fasta_path.empty()) throw UsageError("build: empty FASTA path");
  if (options.chunk_bits == 0) throw UsageError("build: chunk_bits must be > 0");
  const unsigned workers = options.threads == 0 ? 1u : options.threads;

  // Source metadata first: fails fast and is the only hashing work.
  std::error_code ec;
  const uint64_t file_bytes = fs::file_size(fasta_path, ec);
  if (ec) throw InputError("cannot stat '" + fasta_path + "': " + ec.message());
  std::string digest;
  try {
    digest = sha256_file_hex(fasta_path);
  } catch (const IoError& error) {
    throw InputError(error.what());
  }

  genome_.clear();
  valid_.clear();
  contigs_.clear();
  manifest_ = Manifest();
  built_ = false;

  const auto started = std::chrono::steady_clock::now();
  if (options.verbose)
    std::print(stderr,
               "[gcrank] building from {} ({} bytes) with {} worker thread(s), "
               "chunk {} bases\n",
               fasta_path, file_bytes, workers, options.chunk_bits);

  BuildPipeline pipeline(fasta_path, options.chunk_bits, workers, options.verbose);
  pipeline.run(genome_, valid_, contigs_);

  genome_.finalize();
  valid_.finalize();

  manifest_.format_version = kFormatVersion;
  manifest_.builder_version = kBuilderVersion;
  manifest_.created = created_from_env();
  manifest_.source.path = fasta_path;
  manifest_.source.sha256 = digest;
  manifest_.source.file_bytes = file_bytes;
  manifest_.source.total_bases = genome_.size();
  manifest_.source.num_contigs = static_cast<uint64_t>(contigs_.size());
  manifest_.bit_order = kBitOrderMsb0;
  manifest_.endianness = host_endianness_name();
  manifest_.superblock_bits = kSuperblockBits;
  manifest_.basic_block_bits = kBasicBlockBits;
  manifest_.remove_ambiguous_supported = true;
  manifest_.semantics = default_semantics();
  built_ = true;

  if (options.verbose) {
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    std::print(stderr, "[gcrank] built {} bases, {} contigs, {} GC in {:.2f}s\n",
               manifest_.source.total_bases, manifest_.source.num_contigs,
               genome_.rank1(genome_.size()), seconds);
  }
}

void GcRank::save(const std::string& out_dir, const std::string& name, bool overwrite) const {
  if (!built_) throw Error("save: no index in memory (call build() or load() first)");
  if (name.empty()) throw UsageError("save: index name must not be empty");
  if (name == "." || name == "..") throw UsageError("save: invalid index name '" + name + "'");
  if (name.find('/') != std::string::npos || name.find('\\') != std::string::npos)
    throw UsageError("save: index name must not contain path separators: '" + name + "'");
  if (!genome_.has_rank() || !valid_.has_rank())
    throw Error("save: rank indices are not built");

  const fs::path out(out_dir);
  const fs::path final_dir = out / (name + ".gcrank");
  const fs::path tmp_dir = out / (name + ".gcrank.tmp");

  std::error_code ec;
  fs::create_directories(out, ec);
  if (ec)
    throw OutputError("cannot create output directory '" + out.string() + "': " + ec.message());
  if (fs::exists(final_dir, ec)) {
    if (!overwrite)
      throw OutputError("'" + final_dir.string() +
                        "' already exists (pass overwrite/--overwrite to replace it)");
    fs::remove_all(final_dir, ec);
    if (ec) throw OutputError("cannot remove '" + final_dir.string() + "': " + ec.message());
  }
  if (fs::exists(tmp_dir, ec)) {
    fs::remove_all(tmp_dir, ec);
    if (ec) throw OutputError("cannot remove leftover '" + tmp_dir.string() + "': " + ec.message());
  }
  fs::create_directories(tmp_dir, ec);
  if (ec) throw OutputError("cannot create '" + tmp_dir.string() + "': " + ec.message());

  try {
    Manifest manifest = manifest_;
    if (manifest.created.empty()) manifest.created = created_from_env();
    manifest.files.clear();
    manifest.files[kContigsFileName] = write_contigs(tmp_dir / kContigsFileName, contigs_);
    manifest.files[kGenomeBvFileName] = genome_.save_bv((tmp_dir / kGenomeBvFileName).string());
    manifest.files[kGenomeRbFileName] = genome_.save_rank((tmp_dir / kGenomeRbFileName).string());
    manifest.files[kValidBvFileName] = valid_.save_bv((tmp_dir / kValidBvFileName).string());
    manifest.files[kValidRbFileName] = valid_.save_rank((tmp_dir / kValidRbFileName).string());
    for (const char* data_file : kDataFiles)
      write_sidecar(tmp_dir / (std::string(data_file) + kCrcSuffix),
                    manifest.files[data_file].crc32);
    manifest.save((tmp_dir / kManifestFileName).string());
    if (const std::string problem = verify_index_files(tmp_dir, manifest); !problem.empty())
      throw OutputError(problem + " after writing; the index was not published");
    fs::rename(tmp_dir, final_dir, ec);
    if (ec)
      throw OutputError("cannot rename '" + tmp_dir.string() + "' to '" + final_dir.string() +
                        "': " + ec.message());
  } catch (...) {
    std::error_code ignored;
    fs::remove_all(tmp_dir, ignored);
    throw;
  }
}

void GcRank::load(const std::string& dir, const OpenOptions& options) {
  const fs::path root(dir);
  Manifest manifest = Manifest::load((root / kManifestFileName).string());
  std::vector<ContigInfo> contigs =
      parse_contigs(root / kContigsFileName, manifest.source.num_contigs,
                    manifest.source.total_bases);
  if (options.validate_crc)
    if (const std::string problem = verify_index_files(root, manifest); !problem.empty())
      throw FormatError(problem);

  const bool stored_little_endian = manifest.endianness == kEndiannessLittle;
  BitVector genome;
  BitVector valid;
  load_bitvector(genome, root / kGenomeBvFileName, root / kGenomeRbFileName,
                 manifest.source.total_bases, stored_little_endian, options.mmap);
  load_bitvector(valid, root / kValidBvFileName, root / kValidRbFileName,
                 manifest.source.total_bases, stored_little_endian, options.mmap);

  genome_ = std::move(genome);
  valid_ = std::move(valid);
  contigs_ = std::move(contigs);
  manifest_ = std::move(manifest);
  built_ = true;
}

void GcRank::check_range(uint64_t start, uint64_t end) const {
  if (!built_) throw Error("no index in memory (call build() or load() first)");
  if (start > end)
    throw QueryError("invalid window [" + std::to_string(start) + ", " + std::to_string(end) +
                     "): start > end");
  if (end > genome_.size())
    throw QueryError("invalid window [" + std::to_string(start) + ", " + std::to_string(end) +
                     "): end exceeds the genome length " + std::to_string(genome_.size()));
}

uint64_t GcRank::gc_count(uint64_t start, uint64_t end) const {
  check_range(start, end);
  return genome_.rank1(end) - genome_.rank1(start);
}

uint64_t GcRank::valid_count(uint64_t start, uint64_t end) const {
  check_range(start, end);
  return valid_.rank1(end) - valid_.rank1(start);
}

double GcRank::gc_fraction(uint64_t start, uint64_t end, bool remove_ambiguous) const {
  check_range(start, end);
  const uint64_t gc = gc_count(start, end);
  const uint64_t length = end - start;
  uint64_t denominator = length;
  if (remove_ambiguous) {
    denominator = valid_count(start, end);
    if (denominator == 0)
      throw QueryError("gc_fraction(remove_ambiguous = true): window [" +
                       std::to_string(start) + ", " + std::to_string(end) +
                       ") holds no unambiguous A/C/G/T bases");
  } else if (length == 0) {
    throw QueryError("gc_fraction: empty window [" + std::to_string(start) + ", " +
                     std::to_string(end) + ")");
  }
  // `bedtools nuc` converts the counts to float and divides in single
  // precision; doing the same makes "%.6f" of the result identical to its
  // pct_gc column.
  const float quotient = static_cast<float>(gc) / static_cast<float>(denominator);
  return static_cast<double>(quotient);
}

uint64_t GcRank::total_bases() const { return genome_.size(); }

const std::vector<ContigInfo>& GcRank::contigs() const { return contigs_; }

const Manifest& GcRank::manifest() const { return manifest_; }

const ContigInfo* GcRank::contig(std::string_view name) const {
  const auto found =
      std::ranges::find_if(contigs_, [name](const ContigInfo& contig) {
        return contig.name == name;
      });
  return found == contigs_.end() ? nullptr : &*found;
}

}  // namespace gcrank
