#include <doctest.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

#include "gcrank/crc32.h"
#include "gcrank/error.h"
#include "gcrank/format.h"
#include "gcrank/gc_rank.h"
#include "gcrank/manifest.h"
#include "gcrank/sha256.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

// The epoch the committed golden index was built with.
constexpr const char* kGoldenEpoch = "1700000000";

void check_against_reference(const GcRank& index, const gctest::NaiveFasta& reference) {
  const uint64_t total = reference.sequence.size();
  CHECK(index.total_bases() == total);
  REQUIRE(index.contigs().size() == reference.names.size());
  for (std::size_t i = 0; i < reference.names.size(); ++i) {
    CHECK(index.contigs()[i].name == reference.names[i]);
    CHECK(index.contigs()[i].length == reference.lengths[i]);
    CHECK(index.contigs()[i].offset == reference.offsets[i]);
  }

  for (const gctest::ExpectedRow& row : gctest::load_expected_windows()) {
    CAPTURE(row.label);
    CHECK(index.gc_count(row.start, row.end) == row.gc_count);
    CHECK(index.valid_count(row.start, row.end) == row.valid_count);
    if (row.gc_fraction == "ERR") {
      CHECK_THROWS_AS(index.gc_fraction(row.start, row.end, false), QueryError);
    } else {
      // Bit-for-bit against the single-precision value implied by the counts.
      CHECK(index.gc_fraction(row.start, row.end, false) ==
            gctest::expect_single(row.gc_count, row.length()));
    }
    if (row.gc_fraction_remove_ambiguous == "ERR") {
      CHECK_THROWS_AS(index.gc_fraction(row.start, row.end, true), QueryError);
    } else {
      CHECK(index.gc_fraction(row.start, row.end, true) ==
            gctest::expect_single(row.gc_count, row.valid_count));
    }
  }

  // Random windows against the naive scan.
  uint64_t state = 0x9E3779B97F4A7C15ull;
  auto draw = [&state] {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return state >> 33;
  };
  for (int i = 0; i < 400; ++i) {
    const uint64_t a = draw() % (total + 1);
    const uint64_t b = draw() % (total + 1);
    const uint64_t start = std::min(a, b);
    const uint64_t end = std::max(a, b);
    CHECK(index.gc_count(start, end) == gctest::naive_gc_count(reference, start, end));
    CHECK(index.valid_count(start, end) == gctest::naive_valid_count(reference, start, end));
    if (end > start) {
      CHECK(index.gc_fraction(start, end, false) ==
            gctest::expect_single(gctest::naive_gc_count(reference, start, end), end - start));
    }
  }
}

}  // namespace

TEST_CASE("roundtrip: build, save and load agree with the reference values") {
  gctest::TempDir tmp("roundtrip");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  const fs::path fasta = gctest::data_file("test_small.fa");
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(fasta);

  GcRank index;
  BuildOptions options;
  options.threads = 2;
  options.chunk_bits = 997;  // many blocks, misaligned with the rank blocks
  index.build(fasta.string(), options);

  check_against_reference(index, reference);
  CHECK(index.contig("chr_n_rich") != nullptr);
  CHECK(index.contig("chr_n_rich")->length == 1000);
  CHECK(index.contig("missing") == nullptr);

  index.save(tmp.path().string(), "test_small");

  GcRank loaded;
  loaded.load(tmp.file("test_small.gcrank").string());
  check_against_reference(loaded, reference);
  CHECK(loaded.manifest().source.sha256 == sha256_file_hex(fasta.string()));
}

TEST_CASE("roundtrip: whole-genome and per-contig totals") {
  GcRank index;
  index.load(gctest::data_file("test_small.gcrank").string());
  for (const gctest::ExpectedTotals& row : gctest::load_expected_totals()) {
    CAPTURE(row.scope);
    uint64_t start = 0;
    uint64_t end = index.total_bases();
    if (row.scope != "genome") {
      REQUIRE(row.scope.rfind("contig:", 0) == 0);
      const std::string name = row.scope.substr(7);
      const ContigInfo* contig = index.contig(name);
      REQUIRE(contig != nullptr);
      start = contig->offset;
      end = contig->offset + contig->length;
    }
    CHECK(end - start == row.length);
    CHECK(index.gc_count(start, end) == row.gc_count);
    CHECK(index.valid_count(start, end) == row.valid_count);
    CHECK(index.gc_fraction(start, end, false) == gctest::expect_single(row.gc_count, row.length));
    CHECK(index.gc_fraction(start, end, true) ==
          gctest::expect_single(row.gc_count, row.valid_count));
  }
}

TEST_CASE("roundtrip: gzip and plain builds produce identical data files") {
  gctest::TempDir tmp("roundtrip-gz");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  GcRank from_plain;
  from_plain.build(gctest::data_file("test_small.fa").string());
  from_plain.save(tmp.path().string(), "plain");

  GcRank from_gzip;
  from_gzip.build(gctest::data_file("test_small.fa.gz").string());
  from_gzip.save(tmp.path().string(), "gzip");

  const fs::path plain_dir = tmp.file("plain.gcrank");
  const fs::path gzip_dir = tmp.file("gzip.gcrank");
  for (const char* name : kDataFiles) {
    CAPTURE(name);
    CHECK(gctest::same_bytes(plain_dir / name, gzip_dir / name));
  }
  // The manifest differs only in the recorded source.
  CHECK(from_plain.manifest().source.sha256 != from_gzip.manifest().source.sha256);
  CHECK(from_plain.total_bases() == from_gzip.total_bases());
  CHECK(from_plain.gc_count(0, from_plain.total_bases()) ==
        from_gzip.gc_count(0, from_gzip.total_bases()));
}

TEST_CASE("roundtrip: a fresh build is byte-identical to the golden index") {
  gctest::TempDir tmp("golden");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  const gctest::CwdGuard cwd(gctest::data_dir());

  GcRank index;
  BuildOptions options;
  options.threads = 1;
  index.build("test_small.fa", options);
  index.save(tmp.path().string(), "test_small");

  const std::string difference =
      gctest::diff_dirs(tmp.file("test_small.gcrank"), gctest::data_file("test_small.gcrank"));
  CHECK(difference.empty());
  CHECK_FALSE(fs::exists(tmp.file("test_small.gcrank.tmp")));
}

TEST_CASE("roundtrip: the golden manifest describes the index completely") {
  const fs::path dir = gctest::data_file("test_small.gcrank");
  GcRank index;
  index.load(dir.string());
  const Manifest& manifest = index.manifest();

  CHECK(manifest.format_version == kFormatVersion);
  CHECK(manifest.builder_version == kBuilderVersion);
  CHECK(manifest.created == "2023-11-14T22:13:20Z");
  CHECK(manifest.source.path == "test_small.fa");
  CHECK(manifest.source.sha256 == sha256_file_hex(gctest::data_file("test_small.fa").string()));
  CHECK(manifest.source.file_bytes == fs::file_size(gctest::data_file("test_small.fa")));
  CHECK(manifest.source.num_contigs == index.contigs().size());
  CHECK(manifest.source.total_bases == index.total_bases());
  CHECK(manifest.bit_order == "msb0");
  CHECK(manifest.endianness == host_endianness_name());
  CHECK(manifest.superblock_bits == 4096);
  CHECK(manifest.basic_block_bits == 512);
  CHECK(manifest.remove_ambiguous_supported);

  CHECK(manifest.files.size() == 5);
  for (const char* name : kDataFiles) {
    CAPTURE(name);
    REQUIRE(manifest.files.count(name) == 1);
    const FileInfo& info = manifest.files.at(name);
    const fs::path path = dir / name;
    CHECK(info.size == fs::file_size(path));
    CHECK(info.crc32 == crc32_file(path.string()));
    // The sidecar repeats the checksum so external tools can verify files
    // without parsing JSON.
    CHECK(gctest::read_text(dir / (std::string(name) + kCrcSuffix)) ==
          crc32_hex(info.crc32) + "\n");
  }

  CHECK(manifest.files.at(kGenomeBvFileName).size == bv_bytes(index.total_bases()));
  CHECK(manifest.files.at(kValidBvFileName).size == bv_bytes(index.total_bases()));
  CHECK(manifest.files.at(kGenomeRbFileName).size == rank_bytes(index.total_bases()));
  CHECK(manifest.files.at(kValidRbFileName).size == rank_bytes(index.total_bases()));

  // contigs.tsv is exactly `<name>\t<length>\t<offset>\n` per contig.
  uint64_t expected = 0;
  for (const ContigInfo& contig : index.contigs())
    expected += contig.name.size() + std::to_string(contig.length).size() +
                std::to_string(contig.offset).size() + 3;
  CHECK(manifest.files.at(kContigsFileName).size == expected);
  CHECK(gctest::read_text(dir / kContigsFileName).rfind('\n') == expected - 1);
}

TEST_CASE("roundtrip: rank overhead stays below the documented 5% per bitvector") {
  // 4096-bit superblocks need 8 bytes and 512-bit basic blocks 2 bytes: three
  // bytes per 512 bitvector bytes (0.585% of the bits), i.e. 4.7% of the
  // bitvector files. Small vectors pay extra for the trailing sentinel entries.
  const uint64_t large = uint64_t(1) << 30;
  const double asymptotic =
      static_cast<double>(rank_bytes(large)) / static_cast<double>(bv_bytes(large));
  CHECK(asymptotic > 0.046);
  CHECK(asymptotic < 0.047);
  const fs::path dir = gctest::data_file("test_small.gcrank");
  const Manifest manifest = Manifest::load((dir / kManifestFileName).string());
  CHECK(fs::file_size(dir / kGenomeRbFileName) == rank_bytes(manifest.source.total_bases));
  CHECK(fs::file_size(dir / kGenomeBvFileName) == bv_bytes(manifest.source.total_bases));
}

TEST_CASE("roundtrip: ambiguous bases are indexed but never counted") {
  gctest::TempDir tmp("iupac");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  const fs::path fasta = tmp.file("iupac.fa");

  // 8 plain bases, then every ambiguity code 10 times in upper and 10 times in
  // lower case (220 bases), then 40 uracils.
  const std::string codes = "NRYSWKMBDHV";
  const uint64_t plain = 8;
  const uint64_t coded = codes.size() * 20;
  const uint64_t uracil = 40;
  std::string sequence = "ACGTACGT";
  for (int round = 0; round < 2; ++round)
    for (int repeat = 0; repeat < 10; ++repeat)
      for (const char c : codes)
        sequence.push_back(round == 0 ? c : static_cast<char>(c - 'A' + 'a'));
  sequence.append(uracil, 'U');
  gctest::write_text(fasta, ">codes\n" + sequence + "\n");

  GcRank index;
  BuildOptions options;
  options.threads = 2;
  options.chunk_bits = 37;  // splits inside the coded run
  index.build(fasta.string(), options);

  const uint64_t total = index.total_bases();
  REQUIRE(total == sequence.size());
  CHECK(total == plain + coded + uracil);
  CHECK(index.contigs().size() == 1);

  // The plain prefix: 4 GC out of 8, whichever denominator is used.
  CHECK(index.gc_count(0, plain) == 4);
  CHECK(index.valid_count(0, plain) == 8);
  CHECK(index.gc_fraction(0, plain, false) == 0.5);
  CHECK(index.gc_fraction(0, plain, true) == 0.5);

  // Codes are neither GC nor valid; uracil is valid (as T) but not GC.
  CHECK(index.gc_count(0, total) == 4);
  CHECK(index.valid_count(0, total) == 8 + uracil);
  CHECK(index.gc_fraction(0, total, false) == gctest::expect_single(4, plain + coded + uracil));
  CHECK(index.gc_fraction(0, total, true) == gctest::expect_single(4, 8 + uracil));

  // A window made only of ambiguity codes. Note that S (G/C), B (C/G/T) and
  // V (A/C/G) are possible GCs and still not counted.
  const uint64_t coded_begin = plain;
  const uint64_t coded_end = plain + coded;
  CHECK(index.gc_count(coded_begin, coded_end) == 0);
  CHECK(index.valid_count(coded_begin, coded_end) == 0);
  CHECK(index.gc_fraction(coded_begin, coded_end, false) == 0.0);
  CHECK_THROWS_AS(index.gc_fraction(coded_begin, coded_end, true), QueryError);

  // A window made only of uracil: valid, never GC, no error either way.
  CHECK(index.gc_count(coded_end, total) == 0);
  CHECK(index.valid_count(coded_end, total) == uracil);
  CHECK(index.gc_fraction(coded_end, total, false) == 0.0);
  CHECK(index.gc_fraction(coded_end, total, true) == 0.0);

  // A window straddling the plain/coded boundary: "ACGT" plus 20 codes.
  const uint64_t begin = 4;
  const uint64_t finish = plain + 20;
  CHECK(finish - begin == 24);
  CHECK(index.gc_count(begin, finish) == 2);
  CHECK(index.valid_count(begin, finish) == 4);
  CHECK(index.gc_fraction(begin, finish, false) == gctest::expect_single(2, 24));
  CHECK(index.gc_fraction(begin, finish, true) == 0.5);

  // The bits survive a save/load round trip and the manifest records that the
  // remove_ambiguous convention is available.
  index.save(tmp.path().string(), "iupac");
  GcRank loaded;
  loaded.load(tmp.file("iupac.gcrank").string());
  CHECK(loaded.manifest().remove_ambiguous_supported);
  CHECK(loaded.total_bases() == total);
  CHECK(loaded.gc_count(0, total) == 4);
  CHECK(loaded.valid_count(0, total) == 8 + uracil);
  CHECK(loaded.gc_count(begin, finish) == 2);
  CHECK_THROWS_AS(loaded.gc_fraction(coded_begin, coded_end, true), QueryError);
  CHECK(loaded.gc_fraction(coded_end, total, true) == 0.0);
}

TEST_CASE("roundtrip: remove_ambiguous is the classic denominator for A/C/G/T/N input") {
  // The original plan's convention: without ambiguity codes in the file, the
  // flag must produce GC / (A+C+G+T) and GC / (end - start) verbatim.
  gctest::TempDir tmp("denominator");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  const fs::path fasta = tmp.file("acgtn.fa");

  // Deterministic 500-base contig with 30 % N.
  std::string sequence;
  uint64_t state = 0x9E3779B97F4A7C15ull;
  for (int i = 0; i < 500; ++i) {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    const char base = "ACGT"[(state >> 33) % 4];
    sequence.push_back((i % 10) < 3 ? 'N' : base);
  }
  gctest::write_text(fasta, ">acgtn\n" + sequence + "\n");

  GcRank index;
  index.build(fasta.string(), BuildOptions{97, 2, false});

  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(fasta);
  const uint64_t total = index.total_bases();
  const uint64_t gc = index.gc_count(0, total);
  const uint64_t acgt = index.valid_count(0, total);
  CHECK(total == 500);
  CHECK(gc == gctest::naive_gc_count(reference, 0, total));
  CHECK(acgt == gctest::naive_valid_count(reference, 0, total));
  CHECK(acgt < total);  // the Ns are not counted

  // The two formulas of the plan, literally, at gcrank's precision.
  CHECK(index.gc_fraction(0, total, false) == gctest::expect_single(gc, total));
  CHECK(index.gc_fraction(0, total, true) == gctest::expect_single(gc, acgt));
  CHECK(acgt == 350);  // 70 % of 500 bases are A/C/G/T

  // And the same holds on the committed fixture contig that carries 30 % N.
  GcRank golden;
  golden.load(gctest::data_file("test_small.gcrank").string());
  const ContigInfo* n_rich = golden.contig("chr_n_rich");
  REQUIRE(n_rich != nullptr);
  const uint64_t begin = n_rich->offset;
  const uint64_t finish = begin + n_rich->length;
  CHECK(golden.valid_count(begin, finish) == 700);
  CHECK(golden.gc_count(begin, finish) == 364);
  CHECK(golden.gc_fraction(begin, finish, false) == gctest::expect_single(364, 1000));
  CHECK(golden.gc_fraction(begin, finish, true) == gctest::expect_single(364, 700));
}

TEST_CASE("roundtrip: gc_fraction is computed like bedtools nuc (single precision)") {
  gctest::TempDir tmp("precision");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  const fs::path fasta = tmp.file("precision.fa");

  // Counts from a real bedtools nuc comparison: this window is one of the 155
  // (of 21000) whose exact quotient rounds differently at six decimals, and
  // `bedtools nuc` prints 0.453301 for it.
  const uint64_t gc = 273711;
  const uint64_t length = 603818;
  std::string sequence(gc, 'G');
  sequence.append(length - gc, 'A');
  gctest::write_text(fasta, ">w\n" + sequence + "\n");

  GcRank index;
  index.build(fasta.string(), BuildOptions{100000, 2, false});
  REQUIRE(index.total_bases() == length);
  REQUIRE(index.gc_count(0, length) == gc);

  const double value = index.gc_fraction(0, length, false);
  CHECK(value == gctest::expect_single(gc, length));
  // ... and not the exact quotient, which would print as 0.453300.
  CHECK(value != static_cast<double>(gc) / static_cast<double>(length));

  char text[32];
  std::snprintf(text, sizeof(text), "%.6f", value);
  CHECK(std::string(text) == "0.453301");

  // The remove_ambiguous denominator is computed the same way (here the window
  // has no ambiguous base, so both denominators are the window length).
  CHECK(index.gc_fraction(0, length, true) == gctest::expect_single(gc, length));
}

TEST_CASE("roundtrip: query arguments are validated") {
  GcRank index;
  index.load(gctest::data_file("test_small.gcrank").string());
  const uint64_t total = index.total_bases();

  CHECK(index.gc_count(3, 3) == 0);
  CHECK(index.valid_count(3, 3) == 0);
  CHECK_THROWS_AS(index.gc_count(10, 5), QueryError);
  CHECK_THROWS_AS(index.gc_count(0, total + 1), QueryError);
  CHECK_THROWS_AS(index.valid_count(0, total + 1), QueryError);
  CHECK_THROWS_AS(index.gc_fraction(5, 5, false), QueryError);  // zero length

  const ContigInfo* n_rich = index.contig("chr_n_rich");
  REQUIRE(n_rich != nullptr);
  // The fixture starts chr_n_rich with three Ns at positions i % 10 < 3.
  CHECK(index.gc_count(n_rich->offset, n_rich->offset + 3) == 0);
  CHECK(index.valid_count(n_rich->offset, n_rich->offset + 3) == 0);
  CHECK(index.gc_fraction(n_rich->offset, n_rich->offset + 3, false) == 0.0);
  CHECK_THROWS_AS(index.gc_fraction(n_rich->offset, n_rich->offset + 3, true), QueryError);

  GcRank empty;
  CHECK_THROWS_AS(empty.gc_count(0, 1), Error);
  CHECK_THROWS_AS(empty.gc_fraction(0, 1), Error);
}

TEST_CASE("roundtrip: corrupted and incomplete indexes are rejected on load") {
  gctest::TempDir tmp("corrupt");
  const fs::path golden = gctest::data_file("test_small.gcrank");

  const auto copy_golden = [&](const std::string& name) {
    const fs::path target = tmp.file(name);
    fs::remove_all(target);
    fs::create_directories(target);
    for (const std::string& file : gctest::list_files(golden)) {
      fs::create_directories((target / file).parent_path());
      fs::copy_file(golden / file, target / file, fs::copy_options::overwrite_existing);
    }
    return target;
  };
  const auto expect_format_error = [](const fs::path& dir, const OpenOptions& options) {
    GcRank index;
    CHECK_THROWS_AS(index.load(dir.string(), options), FormatError);
  };

  {  // A flipped byte in the bitvector is caught by the checksum.
    const fs::path dir = copy_golden("flipped");
    std::vector<uint8_t> bytes = gctest::read_bytes(dir / kGenomeBvFileName);
    REQUIRE(!bytes.empty());
    bytes[0] ^= 0xFF;
    gctest::write_text(dir / kGenomeBvFileName, std::string(bytes.begin(), bytes.end()));
    expect_format_error(dir, OpenOptions{true, false});
    // ... but an unchecked load accepts it.
    GcRank unchecked;
    CHECK_NOTHROW(unchecked.load(dir.string(), OpenOptions{false, false}));
  }
  {  // A truncated contigs.tsv.
    const fs::path dir = copy_golden("short_contigs");
    const std::string text = gctest::read_text(dir / kContigsFileName);
    gctest::write_text(dir / kContigsFileName, text.substr(0, text.size() / 2));
    expect_format_error(dir, OpenOptions{true, false});
  }
  {  // A missing data file.
    const fs::path dir = copy_golden("missing_valid");
    fs::remove(dir / kValidRbFileName);
    expect_format_error(dir, OpenOptions{true, false});
  }
  {  // A missing manifest.
    GcRank index;
    CHECK_THROWS_AS(index.load(tmp.file("nothing-here").string()), IoError);
  }
  {  // A future format version.
    const fs::path dir = copy_golden("future");
    std::string text = gctest::read_text(dir / kManifestFileName);
    const std::string needle = "\"format_version\": \"" + std::string(kFormatVersion) + "\"";
    const std::size_t position = text.find(needle);
    REQUIRE(position != std::string::npos);
    text.replace(position, needle.size(), "\"format_version\": \"2.0.0\"");
    gctest::write_text(dir / kManifestFileName, text);
    expect_format_error(dir, OpenOptions{true, false});
  }
  {  // contigs.tsv whose lengths no longer add up.
    const fs::path dir = copy_golden("bad_offsets");
    std::string text = gctest::read_text(dir / kContigsFileName);
    const std::size_t position = text.find("\t100\t0\n");
    REQUIRE(position != std::string::npos);
    text.replace(position, 8, "\t99\t0\n");
    gctest::write_text(dir / kContigsFileName, text);
    expect_format_error(dir, OpenOptions{true, false});
  }
}

TEST_CASE("roundtrip: build and save failures are typed") {
  gctest::TempDir tmp("failures");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  GcRank index;
  CHECK_THROWS_AS(index.build(tmp.file("absent.fa").string()), InputError);
  CHECK_THROWS_AS(index.build(gctest::data_file("test_small.fa").string(),
                              BuildOptions{0, 1, false}),
                  UsageError);

  const fs::path broken = tmp.file("broken.fa");
  gctest::write_text(broken, ">a\nACGT?\n");
  CHECK_THROWS_AS(index.build(broken.string()), InputError);
  CHECK_THROWS_AS(index.save(tmp.path().string(), "nothing"), Error);  // nothing built

  index.build(gctest::data_file("test_small.fa").string());
  CHECK_THROWS_AS(index.save(tmp.path().string(), ""), UsageError);
  CHECK_THROWS_AS(index.save(tmp.path().string(), "../escape"), UsageError);
  CHECK_THROWS_AS(index.save(tmp.path().string(), "a/b"), UsageError);

  index.save(tmp.path().string(), "small");
  CHECK(fs::exists(tmp.file("small.gcrank")));
  CHECK_THROWS_AS(index.save(tmp.path().string(), "small"), OutputError);
  CHECK_NOTHROW(index.save(tmp.path().string(), "small", true));
  CHECK_FALSE(fs::exists(tmp.file("small.gcrank.tmp")));

  // A save into a path that is a file, not a directory.
  const fs::path blocker = tmp.file("blocker");
  gctest::write_text(blocker, "x");
  CHECK_THROWS_AS(index.save((blocker / "out").string(), "small"), OutputError);
}
