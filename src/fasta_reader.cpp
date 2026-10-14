#include "gcrank/fasta_reader.h"

#include <algorithm>
#include <format>

#include "gcrank/error.h"
#include "gcrank/format.h"

namespace gcrank {

namespace {

// Printable rendering of an offending byte for error messages.
std::string describe_byte(char c) {
  const unsigned char u = static_cast<unsigned char>(c);
  if (u >= 0x20 && u < 0x7f) return std::format("'{}'", c);
  return std::format("0x{:02x}", u);
}

// Extra explanation for the bytes that most often show up in FASTA files that
// are not nucleotide sequences.
std::string invalid_byte_hint(char c) {
  if (c == '-' || c == '.' || c == '*')
    return " (alignment gaps are not nucleotide codes; strip them before indexing)";
  return std::string();
}

}  // namespace

FastaReader::FastaReader(std::unique_ptr<InputStream> in, uint64_t block_bases)
    : in_(std::move(in)), block_bases_(block_bases) {
  if (in_ == nullptr) throw UsageError("FastaReader: null input stream");
  if (block_bases_ == 0) throw UsageError("FastaReader: block size must be > 0");
  lines_ = std::make_unique<LineReader>(*in_);
  // Reserve a sane amount up front; the buffer grows as needed otherwise.
  block_.reserve(static_cast<std::size_t>(std::min<uint64_t>(block_bases_, uint64_t(1) << 20)));
}

FastaReader::~FastaReader() = default;

std::unique_ptr<FastaReader> FastaReader::from_path(const std::string& path,
                                                    uint64_t block_bases) {
  return std::make_unique<FastaReader>(open_input(path), block_bases);
}

bool FastaReader::line_is_blank() const {
  return std::ranges::all_of(
      line_, [](char c) { return base_class(c) == BaseClass::Whitespace; });
}

void FastaReader::parse_header() {
  // The name is the first whitespace-delimited token (samtools faidx style).
  const std::size_t begin = line_.find_first_not_of(kWhitespaceBytes, 1);
  if (begin == std::string_view::npos) {
    next_name_.clear();
    return;
  }
  const std::size_t end = line_.find_first_of(kWhitespaceBytes, begin);
  next_name_.assign(line_, begin,
                    end == std::string_view::npos ? std::string_view::npos : end - begin);
}

FastaReader::AppendResult FastaReader::append_sequence(std::string_view line,
                                                       std::size_t begin) {
  AppendResult result;
  std::size_t i = begin;
  while (i < line.size()) {
    const char c = line[i];
    const BaseClass cls = base_class(c);
    if (cls == BaseClass::Invalid)
      throw InputError(in_->name() + ": invalid character " + describe_byte(c) +
                       " at line " + std::to_string(lines_->line_number()) + ", column " +
                       std::to_string(i + 1) + invalid_byte_hint(c));
    ++i;
    if (cls == BaseClass::Whitespace) continue;
    block_.push_back(normalized_base(c));
    ++contig_bases_;
    ++total_bases_;
    if (block_.size() >= block_bases_) {
      result.consumed = i;
      result.full = true;
      return result;
    }
  }
  result.consumed = i;
  return result;
}

FastaReader::Event FastaReader::next() {
  if (block_valid_) {
    // The caller has consumed the previous block.
    block_.clear();
    block_valid_ = false;
  }

  for (;;) {
    // A finished contig is always reported before anything else.
    if (pending_contig_end_) {
      if (!block_.empty()) {
        block_valid_ = true;
        return Event::Block;  // flush the tail of the contig first
      }
      pending_contig_end_ = false;
      if (at_eof_) phase_ = Phase::Done;
      return Event::ContigEnd;
    }

    // A parsed header waiting to be announced.
    if (header_pending_) {
      header_pending_ = false;
      contig_name_ = next_name_;
      contig_bases_ = 0;
      phase_ = Phase::InSequence;
      ++num_contigs_;
      return Event::ContigBegin;
    }

    // Continue with the unread rest of the current line.
    if (phase_ == Phase::InSequence && have_line_ && line_pos_ < line_.size()) {
      const AppendResult result = append_sequence(line_, line_pos_);
      line_pos_ = result.consumed;
      if (result.full) {
        block_valid_ = true;
        return Event::Block;
      }
      have_line_ = false;
      continue;
    }

    // Read the next line.
    if (!have_line_) {
      const char* data = nullptr;
      std::size_t len = 0;
      if (!lines_->next(&data, &len)) {
        if (phase_ == Phase::InSequence) {
          pending_contig_end_ = true;
          at_eof_ = true;
          continue;
        }
        phase_ = Phase::Done;
        return Event::End;
      }
      line_ = data != nullptr ? std::string_view(data, len) : std::string_view{};
      have_line_ = true;
      line_pos_ = 0;
    }

    // Header line?
    if (!line_.empty() && line_.front() == '>') {
      parse_header();
      have_line_ = false;
      header_pending_ = true;
      if (phase_ == Phase::InSequence) pending_contig_end_ = true;
      continue;
    }

    // Sequence (or blank) line.
    if (phase_ != Phase::InSequence) {
      if (line_is_blank()) {
        have_line_ = false;
        continue;
      }
      throw InputError(in_->name() + ": sequence data at line " +
                       std::to_string(lines_->line_number()) +
                       " before the first '>' header");
    }
    const AppendResult result = append_sequence(line_, line_pos_);
    line_pos_ = result.consumed;
    if (result.full) {
      block_valid_ = true;
      return Event::Block;
    }
    have_line_ = false;
  }
}

}  // namespace gcrank
