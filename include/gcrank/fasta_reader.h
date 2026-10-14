#pragma once
//
// Streaming FASTA parser.
//
// Vocabulary: A/C/G/T (unambiguous DNA), N and the IUPAC ambiguity codes
// R/Y/S/W/K/M/B/D/H/V, plus U (which is folded onto T, as `bedtools nuc` does),
// in either case; the parser normalizes case on the way out. Ambiguous bytes are
// handed to the caller like any other base (they are part of the sequence and of
// the contig length) but they set neither the GC nor the valid bit, so they are
// never counted as G/C or A/C/G/T. Whitespace inside sequence lines is ignored.
// Any other byte (including gaps such as '-', '.' and '*', and a '>' that does
// not start a line) is an InputError. Header lines start with '>' at the first
// byte of a line; the contig name is the first whitespace-delimited token of the
// header, as in `samtools faidx` (the rest of the header line is a description
// and is not stored). Empty contigs and an empty file are accepted.
//
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "gcrank/input_stream.h"
#include "gcrank/line_reader.h"

namespace gcrank {

class FastaReader {
public:
  enum class Event { ContigBegin, Block, ContigEnd, End };

  static constexpr uint64_t kDefaultBlockBases = 100000000ull;

  // Reads from `in`; sequence blocks are at most `block_bases` bases long and
  // never span contigs. Throws UsageError for a null stream or a zero block
  // size, InputError for malformed input.
  FastaReader(std::unique_ptr<InputStream> in,
              uint64_t block_bases = kDefaultBlockBases);
  ~FastaReader();

  // Convenience: opens `path` (plain or gzip, sniffed by magic bytes).
  static std::unique_ptr<FastaReader> from_path(
      const std::string& path, uint64_t block_bases = kDefaultBlockBases);

  // Advances to the next event. The accessors below stay valid until the next
  // call:
  //   ContigBegin -> contig_name() is the name of the new contig
  //   Block       -> block() holds the next <= block_bases bases
  //   ContigEnd   -> contig_bases() is the length of the finished contig
  //   End         -> input exhausted
  Event next();

  const std::string& contig_name() const { return contig_name_; }
  const std::vector<char>& block() const { return block_; }
  uint64_t contig_bases() const { return contig_bases_; }
  uint64_t total_bases() const { return total_bases_; }
  uint64_t num_contigs() const { return num_contigs_; }
  const std::string& source_name() const { return in_->name(); }

private:
  enum class Phase { Start, InSequence, Done };

  // Result of consuming sequence bytes from one line: how far the parser got
  // and whether the block reached its capacity (needs flushing).
  struct AppendResult {
    std::size_t consumed = 0;
    bool full = false;
  };

  void parse_header();
  // Appends sequence bytes of `line` starting at `begin` to the current block.
  AppendResult append_sequence(std::string_view line, std::size_t begin);
  bool line_is_blank() const;

  std::unique_ptr<InputStream> in_;
  std::unique_ptr<LineReader> lines_;
  uint64_t block_bases_;

  Phase phase_ = Phase::Start;
  bool pending_contig_end_ = false;
  bool header_pending_ = false;
  bool at_eof_ = false;
  bool block_valid_ = false;

  std::string_view line_;
  std::size_t line_pos_ = 0;
  bool have_line_ = false;

  std::string next_name_;
  std::string contig_name_;
  std::vector<char> block_;
  uint64_t contig_bases_ = 0;
  uint64_t total_bases_ = 0;
  uint64_t num_contigs_ = 0;
};

}  // namespace gcrank
