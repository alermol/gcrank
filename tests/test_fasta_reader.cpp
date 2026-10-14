#include <doctest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "gcrank/error.h"
#include "gcrank/fasta_reader.h"
#include "gcrank/format.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

struct Collected {
  std::vector<std::string> contig_names;
  std::vector<uint64_t> contig_lengths;
  std::vector<uint64_t> block_sizes;
  std::string sequence;
  size_t contig_events = 0;
  bool saw_end = false;
};

Collected collect(const std::string& path, uint64_t block_bases) {
  Collected out;
  std::unique_ptr<FastaReader> reader = FastaReader::from_path(path, block_bases);
  for (;;) {
    switch (reader->next()) {
      case FastaReader::Event::ContigBegin:
        out.contig_names.push_back(reader->contig_name());
        ++out.contig_events;
        break;
      case FastaReader::Event::Block:
        out.block_sizes.push_back(reader->block().size());
        out.sequence.append(reader->block().begin(), reader->block().end());
        break;
      case FastaReader::Event::ContigEnd:
        out.contig_lengths.push_back(reader->contig_bases());
        ++out.contig_events;
        break;
      case FastaReader::Event::End:
        out.saw_end = true;
        return out;
    }
  }
}

Collected collect(const std::string& path) { return collect(path, 1u << 20); }

void check_block_sizes(const Collected& out, uint64_t block_bases) {
  for (const uint64_t size : out.block_sizes) {
    CHECK(size > 0);
    CHECK(size <= block_bases);
  }
}

}  // namespace

TEST_CASE("fasta: fixture parses the same at every block size") {
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(gctest::data_file("test_small.fa"));
  const uint64_t block_sizes[] = {1, 7, 64, 100, 512, 513, 1024, 4096, 100000000};

  for (const uint64_t block_bases : block_sizes) {
    CAPTURE(block_bases);
    const Collected out = collect(gctest::data_file("test_small.fa").string(), block_bases);
    CHECK(out.saw_end);
    CHECK(out.contig_names == reference.names);
    CHECK(out.contig_lengths == reference.lengths);
    CHECK(out.sequence == reference.sequence);
    CHECK(out.contig_events == reference.names.size() * 2);
    check_block_sizes(out, block_bases);
  }
}

TEST_CASE("fasta: blocks never span contigs") {
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(gctest::data_file("test_small.fa"));
  // Force several blocks per contig.
  Collected out;
  std::unique_ptr<FastaReader> reader =
      FastaReader::from_path(gctest::data_file("test_small.fa").string(), 64);
  uint64_t contig_index = 0;
  uint64_t in_contig = 0;
  for (;;) {
    const FastaReader::Event event = reader->next();
    if (event == FastaReader::Event::End) break;
    if (event == FastaReader::Event::ContigEnd) {
      CHECK(in_contig == reference.lengths[contig_index]);
      ++contig_index;
      in_contig = 0;
    } else if (event == FastaReader::Event::Block) {
      in_contig += reader->block().size();
      CHECK(in_contig <= reference.lengths[contig_index]);
    }
  }
  CHECK(contig_index == reference.names.size());
}

TEST_CASE("fasta: a block size equal to a contig length is respected") {
  gctest::TempDir tmp("fasta-blocksize");
  const fs::path path = tmp.file("small.fa");
  gctest::write_text(path, ">a\nACGTACGTAC\n>b\nACGTAC\n");

  const Collected out = collect(path.string(), 10);
  CHECK(out.block_sizes == std::vector<uint64_t>{10, 6});
  CHECK(out.contig_lengths == std::vector<uint64_t>{10, 6});

  const Collected one_more = collect(path.string(), 11);
  CHECK(one_more.block_sizes == std::vector<uint64_t>{10, 6});
}

TEST_CASE("fasta: CRLF input matches LF input") {
  const Collected lf = collect(gctest::data_file("test_small.fa").string(), 97);
  const Collected crlf = collect(gctest::data_file("test_small_crlf.fa").string(), 97);
  CHECK(crlf.sequence == lf.sequence);
  CHECK(crlf.contig_names == lf.contig_names);
  CHECK(crlf.contig_lengths == lf.contig_lengths);
  CHECK(crlf.block_sizes == lf.block_sizes);
}

TEST_CASE("fasta: gzip is detected by magic bytes, not by extension") {
  gctest::TempDir tmp("fasta-sniff");
  const Collected plain = collect(gctest::data_file("test_small.fa").string(), 128);

  // A gzip stream under a name that does not end in .gz.
  const fs::path disguised = tmp.file("mystery.dat");
  gctest::write_text(disguised, gctest::read_text(gctest::data_file("test_small.fa.gz")));
  const Collected from_gzip = collect(disguised.string(), 128);
  CHECK(from_gzip.sequence == plain.sequence);
  CHECK(from_gzip.contig_names == plain.contig_names);

  // A plain file with a .gz name.
  const fs::path mislabelled = tmp.file("mislabelled.fa.gz");
  gctest::write_text(mislabelled, gctest::read_text(gctest::data_file("test_small.fa")));
  const Collected from_plain = collect(mislabelled.string(), 128);
  CHECK(from_plain.sequence == plain.sequence);
}

TEST_CASE("fasta: empty file and header-only file") {
  gctest::TempDir tmp("fasta-empty");
  const fs::path empty = tmp.file("empty.fa");
  gctest::write_text(empty, "");
  const Collected nothing = collect(empty.string());
  CHECK(nothing.saw_end);
  CHECK(nothing.contig_names.empty());
  CHECK(nothing.block_sizes.empty());

  const fs::path header_only = tmp.file("header_only.fa");
  gctest::write_text(header_only, ">lonely\n");
  const Collected out = collect(header_only.string());
  CHECK(out.saw_end);
  CHECK(out.contig_names == std::vector<std::string>{"lonely"});
  CHECK(out.contig_lengths == std::vector<uint64_t>{0});
  CHECK(out.block_sizes.empty());
}

TEST_CASE("fasta: empty contigs are reported") {
  gctest::TempDir tmp("fasta-emptycontig");
  const fs::path path = tmp.file("gaps.fa");
  gctest::write_text(path, ">a\n>b\nACGT\n>c\n");
  const Collected out = collect(path.string());
  CHECK(out.contig_names == std::vector<std::string>{"a", "b", "c"});
  CHECK(out.contig_lengths == std::vector<uint64_t>{0, 4, 0});
  CHECK(out.sequence == "ACGT");
}

TEST_CASE("fasta: blank lines and embedded whitespace are skipped") {
  gctest::TempDir tmp("fasta-whitespace");
  const fs::path path = tmp.file("spaced.fa");
  gctest::write_text(path,
                     "\n"
                     "   \n"
                     ">chr1 with a description\n"
                     "\n"
                     "AC GT\n"
                     "  ac\tgt  \n"
                     "\n"
                     ">chr2\n"
                     "nnnn\n");
  const Collected out = collect(path.string(), 3);
  CHECK(out.contig_names == std::vector<std::string>{"chr1", "chr2"});
  CHECK(out.contig_lengths == std::vector<uint64_t>{8, 4});
  CHECK(out.sequence == "ACGTACGTNNNN");
  CHECK(out.block_sizes == std::vector<uint64_t>{3, 3, 2, 3, 1});
}

TEST_CASE("fasta: lower case is normalized to upper case") {
  gctest::TempDir tmp("fasta-case");
  const fs::path path = tmp.file("case.fa");
  gctest::write_text(path, ">a\nacgtnACGTN\n");
  const Collected out = collect(path.string());
  CHECK(out.sequence == "ACGTNACGTN");
}

TEST_CASE("format: base classification") {
  // Unambiguous DNA.
  CHECK(base_class('A') == BaseClass::A);
  CHECK(base_class('c') == BaseClass::C);
  CHECK(base_class('G') == BaseClass::G);
  CHECK(base_class('t') == BaseClass::T);
  CHECK(base_is_gc(base_class('G')));
  CHECK(base_is_gc(base_class('c')));
  CHECK(base_is_valid(base_class('T')));
  CHECK(!base_is_valid(base_class('N')));

  // Uracil is the RNA spelling of T: valid, never GC.
  CHECK(base_class('U') == BaseClass::T);
  CHECK(base_class('u') == BaseClass::T);
  CHECK(base_is_valid(base_class('U')));
  CHECK(!base_is_gc(base_class('U')));
  CHECK(normalized_base('u') == 'T');

  // N and every IUPAC ambiguity code: neither GC nor valid.
  const std::string ambiguous = "NRYSWKMBDHVnryswkmbdhv";
  for (const char c : ambiguous) {
    CAPTURE(c);
    CHECK(base_class(c) == BaseClass::Ambiguous);
    CHECK(base_is_ambiguous(base_class(c)));
    CHECK(!base_is_gc(base_class(c)));
    CHECK(!base_is_valid(base_class(c)));
    const char normalized = normalized_base(c);
    const bool upper_cased =
        c < 'a' ? normalized == c : normalized == static_cast<char>(c - 'a' + 'A');
    CHECK(upper_cased);
  }
  // A code that may be G/C is still not counted as GC.
  CHECK(!base_is_gc(base_class('S')));
  CHECK(!base_is_gc(base_class('B')));
  CHECK(!base_is_gc(base_class('V')));

  // Anything else is rejected, whitespace is skipped.
  for (const char c : std::string("XZ*-.0123#")) {
    CAPTURE(c);
    CHECK(base_class(c) == BaseClass::Invalid);
  }

  for (const char c : std::string(" \t\v\f\r")) {
    CAPTURE(static_cast<int>(c));
    CHECK(base_class(c) == BaseClass::Whitespace);
  }
}

TEST_CASE("fasta: IUPAC codes are accepted, U is folded onto T") {
  gctest::TempDir tmp("fasta-iupac");
  const fs::path path = tmp.file("iupac.fa");
  gctest::write_text(path, ">a\n" "ACGTacgt" "nryswkmbdhv" "NRYSWKMBDHV" "Uu" "\n");
  const Collected out = collect(path.string(), 16);
  CHECK(out.sequence == "ACGTACGT" "NRYSWKMBDHV" "NRYSWKMBDHV" "TT");
  CHECK(out.contig_lengths == std::vector<uint64_t>{32});
  CHECK(out.contig_names == std::vector<std::string>{"a"});
  // Codes are ordinary sequence bytes as far as block splitting is concerned.
  CHECK(out.block_sizes == std::vector<uint64_t>{16, 16});
}

TEST_CASE("fasta: very long lines span reader buffers") {
  gctest::TempDir tmp("fasta-longline");
  const fs::path path = tmp.file("long.fa");
  // One sequence line of 200000 bases (well past the 64 KiB read buffer) and a
  // header line that is itself longer than that buffer.
  std::string text = ">long " + std::string(200000, 'x') + "\n";
  for (int i = 0; i < 200000; ++i) text += "ACGT"[i % 4];
  text += "\n";
  gctest::write_text(path, text);

  const Collected out = collect(path.string(), 1u << 16);
  CHECK(out.contig_names == std::vector<std::string>{"long"});
  CHECK(out.contig_lengths == std::vector<uint64_t>{200000});
  CHECK(out.sequence.size() == 200000);
  CHECK(out.sequence.substr(0, 4) == "ACGT");
  CHECK(out.sequence.substr(199996) == "ACGT");
}

TEST_CASE("fasta: no trailing newline") {
  gctest::TempDir tmp("fasta-notrailing");
  const fs::path path = tmp.file("no_eol.fa");
  gctest::write_text(path, ">a\nACGT");
  const Collected out = collect(path.string());
  CHECK(out.sequence == "ACGT");
  CHECK(out.contig_lengths == std::vector<uint64_t>{4});
}

TEST_CASE("fasta: invalid content is rejected") {
  gctest::TempDir tmp("fasta-invalid");

  struct Case {
    const char* name;
    const char* content;
  };
  const Case cases[] = {
      {"invalid_letter.fa", ">a\nACGTXACGT\n"},
      {"digit.fa", ">a\nACGT1\n"},
      {"gt_in_line.fa", ">a\nACGT>ACGT\n"},
      {"sequence_first.fa", "ACGT\n>a\nACGT\n"},
      {"high_byte.fa", ">a\nACGT\xff\n"},
      {"nul.fa", ">a\nACGT\x00ACGT\n"},
      {"gap.fa", ">a\nACGT-ACGT\n"},
      {"stop.fa", ">a\nACGT*ACGT\n"},
  };
  for (const Case& test_case : cases) {
    CAPTURE(test_case.name);
    const fs::path path = tmp.file(test_case.name);
    gctest::write_text(path, test_case.content);
    CHECK_THROWS_AS(collect(path.string()), InputError);
  }
}

TEST_CASE("fasta: missing file and zero block size") {
  CHECK_THROWS_AS(collect("does-not-exist.fa"), InputError);

  gctest::TempDir tmp("fasta-zeroblock");
  const fs::path path = tmp.file("small.fa");
  gctest::write_text(path, ">a\nACGT\n");
  CHECK_THROWS_AS(collect(path.string(), 0), UsageError);
  CHECK_THROWS_AS(FastaReader(nullptr, 10), UsageError);
}

TEST_CASE("fasta: statistics accessors") {
  std::unique_ptr<FastaReader> reader =
      FastaReader::from_path(gctest::data_file("test_small.fa").string(), 512);
  uint64_t contigs = 0;
  uint64_t bases = 0;
  for (;;) {
    const FastaReader::Event event = reader->next();
    if (event == FastaReader::Event::End) break;
    if (event == FastaReader::Event::Block) bases += reader->block().size();
    if (event == FastaReader::Event::ContigEnd) ++contigs;
  }
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(gctest::data_file("test_small.fa"));
  CHECK(contigs == reference.names.size());
  CHECK(bases == reference.sequence.size());
  CHECK(reader->total_bases() == reference.sequence.size());
  CHECK(reader->num_contigs() == reference.names.size());
  CHECK(reader->source_name() == gctest::data_file("test_small.fa").string());
}
