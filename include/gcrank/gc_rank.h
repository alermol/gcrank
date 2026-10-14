#pragma once
//
// GcRank: build, store and query a GC-content index over a FASTA file.
//
// Queries are O(1) in the window length: they are answered from rank indices
// over two bitvectors (GC and valid bases), so a window can cover a whole
// chromosome without scanning it.
//
#include <cstdint>
#include <string>
#include <vector>

#include "gcrank/bitvector.h"
#include "gcrank/manifest.h"

namespace gcrank {

// One contig of the indexed genome.
struct ContigInfo {
  std::string name;
  uint64_t length = 0;  // bases, including N
  uint64_t offset = 0;  // 0-based offset of the first base
};

struct BuildOptions {
  // Bases per parallel work unit. Blocks never span contigs and the number of
  // blocks only affects scheduling, never the bytes written. Each in-flight
  // block costs roughly 1.25x chunk_bits bytes (raw bases plus the two packed
  // bitvector slices), across threads + 2 slots, so lower this for large
  // thread counts on memory-constrained machines.
  uint64_t chunk_bits = 100000000ull;
  unsigned threads = 4;  // worker threads (0 -> 1); the producer and the
                         // assembling caller are extra
  bool verbose = false;  // progress on stderr
};

struct OpenOptions {
  bool validate_crc = true;  // re-read and verify every data file
  bool mmap = false;         // map data files instead of reading them
};

// Default worker count: min(8, hardware_concurrency()), never below 1.
unsigned default_threads();

class GcRank {
public:
  GcRank();
  ~GcRank();
  GcRank(GcRank&&) noexcept;
  GcRank& operator=(GcRank&&) noexcept;
  GcRank(const GcRank&) = delete;
  GcRank& operator=(const GcRank&) = delete;

  // Parses `fasta_path` (plain or gzip) and builds both bitvectors and their
  // rank indices in memory. Parallel: one producer thread, `threads` workers,
  // and the calling thread assembles. The result is byte-identical for any
  // thread count. Throws InputError for unreadable/malformed input.
  void build(const std::string& fasta_path, const BuildOptions& options = {});

  // Writes the index to `<out_dir>/<name>.gcrank/`. The index is written to a
  // temporary directory, verified against its checksums, and then renamed into
  // place, so an interrupted save never leaves a half-written index behind.
  // Refuses to replace an existing index unless `overwrite` is true.
  void save(const std::string& out_dir, const std::string& name,
            bool overwrite = false) const;

  // Loads an index directory (the `<name>.gcrank` directory itself).
  void load(const std::string& dir, const OpenOptions& options = {});

  // ---- queries: half-open [start, end), 0-based --------------------------

  // Number of unambiguous G/C bases in the window. Ambiguity codes (N, U, R,
  // Y, S, W, K, M, B, D, H, V) are never counted, so this is a lower bound for
  // windows that contain them.
  uint64_t gc_count(uint64_t start, uint64_t end) const;
  // Number of unambiguous A/C/G/T bases in the window; N, U and the IUPAC
  // ambiguity codes are excluded.
  uint64_t valid_count(uint64_t start, uint64_t end) const;
  // GC fraction of the window, computed exactly like `bedtools nuc`: the
  // integer counts are converted to `float` and divided in single precision,
  // and the result is returned as a double holding that single-precision
  // value. Printing it with "%.6f" is therefore identical to bedtools' pct_gc
  // column, down to the last digit (above 2^24 bases per window the float
  // conversion rounds the counts themselves, the same way bedtools does).
  //
  // With remove_ambiguous = false the denominator is the window length, i.e.
  // everything the window spans (bedtools' `seq_len`); with
  // remove_ambiguous = true it is valid_count, so every ambiguous base is
  // removed (bedtools has no such column). The window must then hold at least
  // one A/C/G/T base, otherwise QueryError is thrown, as it is for a
  // zero-length window either way.
  double gc_fraction(uint64_t start, uint64_t end, bool remove_ambiguous = false) const;

  // ---- metadata ----------------------------------------------------------

  uint64_t total_bases() const;
  const std::vector<ContigInfo>& contigs() const;
  const Manifest& manifest() const;
  // Linear lookup by contig name; returns nullptr when absent.
  const ContigInfo* contig(const std::string& name) const;

private:
  void check_range(uint64_t start, uint64_t end) const;

  BitVector genome_;
  BitVector valid_;
  std::vector<ContigInfo> contigs_;
  Manifest manifest_;
  bool built_ = false;
};

}  // namespace gcrank
