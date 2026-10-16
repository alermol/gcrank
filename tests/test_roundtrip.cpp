#include <doctest.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#endif

#include "gcrank/detail/crc32.h"
#include "gcrank/error.h"
#include "gcrank/format.h"
#include "gcrank/gc_rank.h"
#include "gcrank/manifest.h"
#include "gcrank/detail/sha256.h"
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
      CHECK_THROWS_AS((void)index.gc_fraction(row.start, row.end, false), QueryError);
    } else {
      // Bit-for-bit against the single-precision value implied by the counts.
      CHECK(index.gc_fraction(row.start, row.end, false) ==
            gctest::expect_single(row.gc_count, row.length()));
    }
    if (row.gc_fraction_remove_ambiguous == "ERR") {
      CHECK_THROWS_AS((void)index.gc_fraction(row.start, row.end, true), QueryError);
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
  // The manifest differs only in the recorded source: the hash covers the
  // parsed byte stream, so a file and its gzip form hash identically (and the
  // recorded size is the decompressed size).
  CHECK(from_plain.manifest().source.sha256 == from_gzip.manifest().source.sha256);
  CHECK(from_plain.manifest().source.file_bytes == from_gzip.manifest().source.file_bytes);
  CHECK(from_gzip.manifest().source.file_bytes ==
        fs::file_size(gctest::data_file("test_small.fa")));
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
  // 1.1 writers store .rb scalars little-endian on every host.
  CHECK(manifest.endianness == kEndiannessLittle);
  CHECK(manifest.has_totals);
  CHECK(manifest.total_gc == index.gc_count(0, index.total_bases()));
  CHECK(manifest.total_valid == index.valid_count(0, index.total_bases()));
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
  CHECK_THROWS_AS((void)index.gc_fraction(coded_begin, coded_end, true), QueryError);

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
  CHECK_THROWS_AS((void)loaded.gc_fraction(coded_begin, coded_end, true), QueryError);
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
  index.build(fasta.string(), BuildOptions{97, 2, false, {}});

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
  index.build(fasta.string(), BuildOptions{100000, 2, false, {}});
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
  CHECK_THROWS_AS((void)index.gc_count(10, 5), QueryError);
  CHECK_THROWS_AS((void)index.gc_count(0, total + 1), QueryError);
  CHECK_THROWS_AS((void)index.valid_count(0, total + 1), QueryError);
  CHECK_THROWS_AS((void)index.gc_fraction(5, 5, false), QueryError);  // zero length

  const ContigInfo* n_rich = index.contig("chr_n_rich");
  REQUIRE(n_rich != nullptr);
  // The fixture starts chr_n_rich with three Ns at positions i % 10 < 3.
  CHECK(index.gc_count(n_rich->offset, n_rich->offset + 3) == 0);
  CHECK(index.valid_count(n_rich->offset, n_rich->offset + 3) == 0);
  CHECK(index.gc_fraction(n_rich->offset, n_rich->offset + 3, false) == 0.0);
  CHECK_THROWS_AS((void)index.gc_fraction(n_rich->offset, n_rich->offset + 3, true), QueryError);

  GcRank empty;
  CHECK_THROWS_AS((void)empty.gc_count(0, 1), Error);
  CHECK_THROWS_AS((void)empty.gc_fraction(0, 1), Error);
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
    expect_format_error(dir, OpenOptions{true, false, {}});
    // ... but an unchecked load accepts it.
    GcRank unchecked;
    CHECK_NOTHROW(unchecked.load(dir.string(), OpenOptions{false, false, {}}));
  }
  {  // A truncated contigs.tsv.
    const fs::path dir = copy_golden("short_contigs");
    const std::string text = gctest::read_text(dir / kContigsFileName);
    gctest::write_text(dir / kContigsFileName, text.substr(0, text.size() / 2));
    expect_format_error(dir, OpenOptions{true, false, {}});
  }
  {  // A missing data file.
    const fs::path dir = copy_golden("missing_valid");
    fs::remove(dir / kValidRbFileName);
    expect_format_error(dir, OpenOptions{true, false, {}});
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
    text.replace(position, needle.size(), "\"format_version\": \"9.0.0\"");
    gctest::write_text(dir / kManifestFileName, text);
    expect_format_error(dir, OpenOptions{true, false, {}});
  }
  {  // contigs.tsv whose lengths do not add up.
    const fs::path dir = copy_golden("bad_offsets");
    std::string text = gctest::read_text(dir / kContigsFileName);
    const std::size_t position = text.find("\t100\t0\n");
    REQUIRE(position != std::string::npos);
    text.replace(position, 8, "\t99\t0\n");
    gctest::write_text(dir / kContigsFileName, text);
    expect_format_error(dir, OpenOptions{true, false, {}});
  }
}

TEST_CASE("roundtrip: build and save failures are typed") {
  gctest::TempDir tmp("failures");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  GcRank index;
  CHECK_THROWS_AS(index.build(tmp.file("absent.fa").string()), InputError);
  CHECK_THROWS_AS(index.build(gctest::data_file("test_small.fa").string(),
                              BuildOptions{0, 1, false, {}}),
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

TEST_CASE("roundtrip: a failed overwrite keeps the old index") {
  gctest::TempDir tmp("overwrite-failure");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  GcRank index;
  index.build(gctest::data_file("test_small.fa").string());
  const fs::path out = tmp.file("out");
  index.save(out.string(), "small");
  const fs::path index_dir = out / "small.gcrank";
  REQUIRE(fs::exists(index_dir / kManifestFileName));
  const std::vector<uint8_t> manifest_before =
      gctest::read_bytes(index_dir / kManifestFileName);
  const std::vector<uint8_t> genome_before = gctest::read_bytes(index_dir / kGenomeBvFileName);

#if defined(__linux__) || defined(__APPLE__)
  // Only meaningful as a non-root user: root bypasses directory permissions.
  if (::geteuid() != 0) {
    fs::permissions(out, fs::perms::owner_read | fs::perms::owner_exec);
    CHECK_THROWS_AS(index.save(out.string(), "small", /*overwrite=*/true), OutputError);
    fs::permissions(out, fs::perms::owner_all, fs::perm_options::replace);
    CHECK(gctest::read_bytes(index_dir / kManifestFileName) == manifest_before);
    CHECK(gctest::read_bytes(index_dir / kGenomeBvFileName) == genome_before);
    GcRank reloaded;
    CHECK_NOTHROW(reloaded.load(index_dir.string()));
    CHECK(reloaded.total_bases() == index.total_bases());
  }
#endif
}

TEST_CASE("roundtrip: concurrent saves of the same name do not corrupt each other") {
  gctest::TempDir tmp("concurrent-save");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  GcRank index;
  BuildOptions options;
  options.threads = 2;
  index.build(gctest::data_file("test_small.fa").string(), options);
  index.save(tmp.path().string(), "ref");  // reference bytes for the same data

  std::exception_ptr first_error;
  std::exception_ptr second_error;
  {
    std::thread a([&] {
      try {
        index.save(tmp.path().string(), "same", /*overwrite=*/true);
      } catch (...) {
        first_error = std::current_exception();
      }
    });
    std::thread b([&] {
      try {
        index.save(tmp.path().string(), "same", /*overwrite=*/true);
      } catch (...) {
        second_error = std::current_exception();
      }
    });
    a.join();
    b.join();
  }
  if (first_error != nullptr || second_error != nullptr) {
    try {
      std::rethrow_exception(first_error != nullptr ? first_error : second_error);
    } catch (const std::exception& error) {
      FAIL("concurrent save threw: " << error.what());
    } catch (...) {
      FAIL("concurrent save threw a non-standard exception");
    }
  }
  CHECK(gctest::diff_dirs(tmp.file("same.gcrank").string(),
                          tmp.file("ref.gcrank").string())
            .empty());
  for (const auto& entry : fs::directory_iterator(tmp.path())) {
    CAPTURE(entry.path().string());
    const std::string name = entry.path().filename().string();
    CHECK(name.find(".tmp.") == std::string::npos);
    CHECK(name.find(".old.") == std::string::npos);
  }
}

TEST_CASE("roundtrip: SOURCE_DATE_EPOCH is validated before the input is read") {
  gctest::TempDir tmp("bad-epoch");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", "not-an-epoch");
  GcRank index;
  // The input does not exist, so an InputError would be the natural failure;
  // the invalid epoch must win so that a multi-hour build cannot fail on it at
  // the very end.
  CHECK_THROWS_AS(index.build(tmp.file("absent.fa").string()), UsageError);
}

TEST_CASE("roundtrip: contig-relative queries match the global ones and reject crossings") {
  gctest::TempDir tmp("contig-relative");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  GcRank index;
  index.build(gctest::data_file("test_small.fa").string());

  for (const gctest::ExpectedTotals& row : gctest::load_expected_totals()) {
    if (!row.scope.starts_with("contig:")) continue;
    const std::string name = row.scope.substr(std::string_view("contig:").size());
    CAPTURE(name);
    const ContigInfo* info = index.contig(name);
    REQUIRE(info != nullptr);
    CHECK(info->length == row.length);
    CHECK(index.gc_count(name, 0, info->length) == row.gc_count);
    CHECK(index.valid_count(name, 0, info->length) == row.valid_count);
    CHECK(index.gc_count(name, 0, info->length) ==
          index.gc_count(info->offset, info->offset + info->length));
    if (info->length >= 4) {
      const uint64_t begin = info->length / 3;
      const uint64_t end = info->length - 2;
      CHECK(index.gc_count(name, begin, end) ==
            index.gc_count(info->offset + begin, info->offset + end));
      CHECK(index.valid_count(name, begin, end) ==
            index.valid_count(info->offset + begin, info->offset + end));
      CHECK(index.gc_fraction(name, begin, end, false) ==
            index.gc_fraction(info->offset + begin, info->offset + end, false));
    }
  }

  // The classic misuse: a BED interval that runs past the contig end must fail
  // instead of silently continuing into the next contig.
  const ContigInfo& first = index.contigs().front();
  CHECK_THROWS_AS(static_cast<void>(index.gc_count(first.name, 0, first.length + 1)), QueryError);
  CHECK_THROWS_AS(static_cast<void>(index.valid_count(first.name, first.length, first.length + 1)),
                  QueryError);
  CHECK_THROWS_AS(static_cast<void>(index.gc_fraction(first.name, 5, 4)), QueryError);
  CHECK_THROWS_AS(static_cast<void>(index.gc_count(std::string_view("not_a_contig"), 0, 1)),
                  QueryError);
  CHECK(index.gc_count(first.name, 0, 0) == 0);
  CHECK_THROWS_AS(static_cast<void>(index.gc_fraction(first.name, 0, 0)), QueryError);

  // contig() is an exact match (no prefix or suffix hits).
  CHECK(index.contig("chr_boundary_100") != nullptr);
  CHECK(index.contig("chr_boundary") == nullptr);
  CHECK(index.contig("chr_boundary_100x") == nullptr);
  CHECK(index.contig("") == nullptr);

  // An empty contig is a valid name with length 0.
  const fs::path empty_contig_fasta = tmp.file("empty_contig.fa");
  gctest::write_text(empty_contig_fasta, ">empty\n>a\nACGT\n");
  GcRank with_empty;
  with_empty.build(empty_contig_fasta.string());
  REQUIRE(with_empty.contig("empty") != nullptr);
  CHECK(with_empty.contig("empty")->length == 0);
  CHECK(with_empty.gc_count("empty", 0, 0) == 0);
  CHECK_THROWS_AS(static_cast<void>(with_empty.gc_count("empty", 0, 1)), QueryError);

  // Queries on an index that was never built or loaded are usage errors.
  GcRank fresh;
  CHECK_THROWS_AS(static_cast<void>(fresh.gc_count(0, 1)), UsageError);
  CHECK_THROWS_AS(static_cast<void>(fresh.gc_count("chr1", 0, 1)), UsageError);
  CHECK(fresh.contig("chr1") == nullptr);
}

TEST_CASE("roundtrip: a failed build leaves the previous index usable") {
  gctest::TempDir tmp("strong-guarantee");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);
  GcRank index;
  index.build(gctest::data_file("test_small.fa").string());
  const uint64_t total = index.total_bases();
  const uint64_t gc = index.gc_count(0, total);
  const std::size_t contigs = index.contigs().size();
  REQUIRE(index.contig("chr_n_rich") != nullptr);

  const fs::path broken = tmp.file("broken.fa");
  gctest::write_text(broken, ">a\nACGT?\n");
  CHECK_THROWS_AS(index.build(broken.string()), InputError);
  CHECK_THROWS_AS(index.build(tmp.file("absent.fa").string()), InputError);

  // The loaded index survived both failures intact.
  CHECK(index.total_bases() == total);
  CHECK(index.gc_count(0, total) == gc);
  CHECK(index.contigs().size() == contigs);
  CHECK(index.contig("chr_n_rich") != nullptr);
  CHECK_NOTHROW(index.save(tmp.path().string(), "still_here"));
}

TEST_CASE("roundtrip: corruption is caught without full checksum validation") {
  gctest::TempDir tmp("structural-checks");
  const fs::path golden = gctest::data_file("test_small.gcrank");
  const auto copy_golden = [&](const std::string& name) {
    const fs::path target = tmp.file(name);
    fs::remove_all(target);
    fs::copy(golden, target, fs::copy_options::recursive);
    return target;
  };
  constexpr uint64_t kGoldenBits = 11395;
  constexpr bool kNoCrc = false;
  constexpr bool kNoMmap = false;

  {  // A tampered sidecar is rejected even when the data CRCs are skipped.
    const fs::path dir = copy_golden("sidecar");
    gctest::write_text(dir / (std::string(kGenomeBvFileName) + kCrcSuffix), "00000000\n");
    GcRank index;
    CHECK_THROWS_AS(index.load(dir.string(), OpenOptions{kNoCrc, kNoMmap, {}}), FormatError);
  }
  {  // A corrupted rank index (a used superblock rank bumped) changes the
     // population total and is rejected through the manifest cross-check.
    const fs::path dir = copy_golden("rank_total");
    std::vector<uint8_t> bytes = gctest::read_bytes(dir / kGenomeRbFileName);
    const uint64_t nsb = rank_superblock_entries(kGoldenBits);
    REQUIRE(bytes.size() == rank_bytes(kGoldenBits));
    // The last used entry is nsb - 2; the trailing one is padding for
    // positions past the end of the vector and is never read by rank1().
    const std::size_t at = static_cast<std::size_t>((nsb - 2) * 8);
    uint64_t value = 0;
    for (int b = 0; b < 8; ++b) value |= uint64_t(bytes[at + static_cast<std::size_t>(b)]) << (8 * b);
    value += 1;
    for (int b = 0; b < 8; ++b)
      bytes[at + static_cast<std::size_t>(b)] = static_cast<uint8_t>((value >> (8 * b)) & 0xFFu);
    gctest::write_text(dir / kGenomeRbFileName, std::string(bytes.begin(), bytes.end()));
    GcRank index;
    CHECK_THROWS_AS(index.load(dir.string(), OpenOptions{kNoCrc, kNoMmap, {}}), FormatError);
  }
  {  // A basic-block entry that does not restart at a superblock boundary is
     // rejected by the rank validator itself.
    const fs::path dir = copy_golden("rank_bb");
    std::vector<uint8_t> bytes = gctest::read_bytes(dir / kGenomeRbFileName);
    const uint64_t nsb = rank_superblock_entries(kGoldenBits);
    const std::size_t at = static_cast<std::size_t>(8 * nsb + 2 * 8);  // bit 4096
    REQUIRE(at + 1 < bytes.size());
    bytes[at] = 1;
    bytes[at + 1] = 0;
    gctest::write_text(dir / kGenomeRbFileName, std::string(bytes.begin(), bytes.end()));
    GcRank index;
    CHECK_THROWS_AS(index.load(dir.string(), OpenOptions{kNoCrc, kNoMmap, {}}), FormatError);
  }
  {  // A manifest without the optional population totals still loads.
    const fs::path dir = copy_golden("no_totals");
    Manifest manifest = Manifest::load((dir / kManifestFileName).string());
    manifest.has_totals = false;
    gctest::write_text(dir / kManifestFileName, manifest.to_json());
    GcRank index;
    CHECK_NOTHROW(index.load(dir.string()));
    CHECK_FALSE(index.manifest().has_totals);
    CHECK(index.gc_count(0, index.total_bases()) > 0);
  }
}

TEST_CASE("roundtrip: the library logs through the caller's sink, never stderr") {
  gctest::TempDir tmp("logging");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", kGoldenEpoch);

  std::mutex mutex;
  std::vector<std::pair<LogLevel, std::string>> messages;
  const LogSink sink = [&](LogLevel level, std::string_view message) {
    std::lock_guard<std::mutex> lock(mutex);
    messages.emplace_back(level, std::string(message));
  };
  const auto take_messages = [&] {
    std::lock_guard<std::mutex> lock(mutex);
    return messages;
  };

  GcRank index;
  BuildOptions options;
  options.verbose = true;
  options.log = sink;
  index.build(gctest::data_file("test_small.fa").string(), options);
  {
    const auto taken = take_messages();
    REQUIRE(!taken.empty());
    CHECK(taken.front().first == LogLevel::Info);
    CHECK(taken.front().second.find("building from") != std::string::npos);
    CHECK(taken.back().first == LogLevel::Info);
    CHECK(taken.back().second.find("built") != std::string::npos);
  }

  // Without verbose the sink only ever sees warnings.
  messages.clear();
  options.verbose = false;
  index.build(gctest::data_file("test_small.fa").string(), options);
  for (const auto& entry : take_messages()) CHECK(entry.first != LogLevel::Info);

  // The mmap fallback is reported as a warning (and the load fails on the
  // missing file afterwards, which is fine: the warning is the subject here).
  const fs::path dir = tmp.file("idx.gcrank");
  index.save(tmp.path().string(), "idx");
  fs::remove(dir / kValidBvFileName);
  messages.clear();
  OpenOptions open;
  open.validate_crc = false;  // reach the mmap path; the CRC pass would fail first
  open.mmap = true;
  open.log = sink;
  GcRank broken;
  CHECK_THROWS_AS(broken.load(dir.string(), open), IoError);
  bool warned = false;
  for (const auto& entry : take_messages()) warned |= entry.first == LogLevel::Warning;
  CHECK(warned);
}
