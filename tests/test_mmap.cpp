#include <doctest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "gcrank/error.h"
#include "gcrank/gc_rank.h"
#include "gcrank/detail/input_stream.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

void check_same_answers(const GcRank& left, const GcRank& right, const fs::path& fasta) {
  REQUIRE(left.total_bases() == right.total_bases());
  const uint64_t total = left.total_bases();
  CHECK(left.gc_count(0, total) == right.gc_count(0, total));
  CHECK(left.valid_count(0, total) == right.valid_count(0, total));

  uint64_t state = 0xC0FFEE123456789ull;
  auto draw = [&state] {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return state >> 33;
  };
  for (int i = 0; i < 200; ++i) {
    const uint64_t a = draw() % (total + 1);
    const uint64_t b = draw() % (total + 1);
    const uint64_t start = std::min(a, b);
    const uint64_t end = std::max(a, b);
    CHECK(left.gc_count(start, end) == right.gc_count(start, end));
    CHECK(left.valid_count(start, end) == right.valid_count(start, end));
    if (start != end) {
      CHECK(left.gc_fraction(start, end, false) == right.gc_fraction(start, end, false));
      if (right.valid_count(start, end) > 0)
        CHECK(left.gc_fraction(start, end, true) == right.gc_fraction(start, end, true));
    }
  }

  // And against the independent scan of the source file.
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(fasta);
  CHECK(left.gc_count(0, total) == gctest::naive_gc_count(reference, 0, total));
  CHECK(left.valid_count(0, total) == gctest::naive_valid_count(reference, 0, total));
  CHECK(left.gc_count(total / 3, total - 7) ==
        gctest::naive_gc_count(reference, total / 3, total - 7));
}

}  // namespace

TEST_CASE("mmap: mapped and heap loads answer identically") {
  gctest::TempDir tmp("mmap-compare");
  const fs::path fasta = tmp.file("synthetic.fa");
  gctest::write_synthetic_fasta(fasta, 4000000, 700000, /*n_every=*/97);

  GcRank builder;
  BuildOptions options;
  options.threads = 4;
  options.chunk_bits = 100000;
  builder.build(fasta.string(), options);
  builder.save(tmp.path().string(), "synthetic");
  const fs::path dir = tmp.file("synthetic.gcrank");

  GcRank heap;
  heap.load(dir.string(), OpenOptions{/*validate_crc=*/true, /*mmap=*/false, {}});
  GcRank mapped;
  mapped.load(dir.string(), OpenOptions{/*validate_crc=*/true, /*mmap=*/true, {}});
  check_same_answers(heap, mapped, fasta);
  check_same_answers(builder, mapped, fasta);

  // Without checksum validation the mapped load must still be usable.
  GcRank mapped_fast;
  mapped_fast.load(dir.string(), OpenOptions{false, true, {}});
  check_same_answers(heap, mapped_fast, fasta);
}

TEST_CASE("mmap: mapped files are not copied into heap memory") {
  if (gctest::rss_bytes() == 0) {
    MESSAGE("RSS measurement is unavailable on this platform; skipping");
    return;
  }
  gctest::TempDir tmp("mmap-rss");
  const fs::path fasta = tmp.file("synthetic.fa");
  const uint64_t total = 24000000;  // ~24 MB, large enough to see page effects
  gctest::write_synthetic_fasta(fasta, total, 2000000);

  const uint64_t bitvector_bytes = 2 * bv_bytes(total);
  {
    GcRank builder;
    BuildOptions options;
    options.threads = 4;
    options.chunk_bits = 4000000;
    builder.build(fasta.string(), options);
    builder.save(tmp.path().string(), "synthetic");
  }
  const fs::path dir = tmp.file("synthetic.gcrank");

  // Parse the reference (and its 24 MB of sequence) before measuring so that
  // it is part of the baseline rather than of a measured delta.
  const gctest::NaiveFasta reference = gctest::parse_fasta_plain(fasta);
  const uint64_t expected_gc = gctest::naive_gc_count(reference, 0, total);
  const uint64_t expected_valid = gctest::naive_valid_count(reference, 0, total);

  gctest::release_free_memory();
  const uint64_t rss_before = gctest::rss_bytes();
  int64_t mmap_growth = 0;
  {
    GcRank mapped;
    mapped.load(dir.string(), OpenOptions{/*validate_crc=*/false, /*mmap=*/true, {}});
    mmap_growth = static_cast<int64_t>(gctest::rss_bytes()) - static_cast<int64_t>(rss_before);
    CHECK(mapped.gc_count(0, total) == expected_gc);
    CHECK(mapped.valid_count(0, total) == expected_valid);
  }

  gctest::release_free_memory();
  const uint64_t rss_after_mmap = gctest::rss_bytes();
  int64_t heap_growth = 0;
  {
    GcRank heap;
    heap.load(dir.string(), OpenOptions{/*validate_crc=*/false, /*mmap=*/false, {}});
    heap_growth = static_cast<int64_t>(gctest::rss_bytes()) - static_cast<int64_t>(rss_after_mmap);
    CHECK(heap.gc_count(0, total) == expected_gc);
    CHECK(heap.valid_count(0, total) == expected_valid);
  }

  MESSAGE("bitvector bytes: " << bitvector_bytes << ", mmap RSS growth: " << mmap_growth
                              << ", heap RSS growth: " << heap_growth);
  // A mapped load only touches the rank arrays (~4.7% of the bitvectors).
  CHECK(mmap_growth < static_cast<int64_t>(bitvector_bytes / 2));
  CHECK(heap_growth > mmap_growth);
}

TEST_CASE("mmap: MappedFile contract") {
  gctest::TempDir tmp("mmap-contract");
  const fs::path path = tmp.file("data.bin");
  const std::string text = "gcrank mapping test";
  gctest::write_text(path, text);

  MappedFile file = MappedFile::map(path.string());
  REQUIRE(file.valid());
  CHECK(file.size() == text.size());
  CHECK(std::string(reinterpret_cast<const char*>(file.data()), file.size()) == text);
  CHECK(file.name() == path.string());

  MappedFile moved = std::move(file);
  CHECK(moved.valid());
  CHECK(moved.size() == text.size());
  moved.reset();
  CHECK_FALSE(moved.valid());
  CHECK(moved.data() == nullptr);

  // Mapping an empty file succeeds with an empty view.
  const fs::path empty = tmp.file("empty.bin");
  gctest::write_text(empty, "");
  MappedFile empty_map = MappedFile::map(empty.string());
  CHECK(empty_map.valid());
  CHECK(empty_map.size() == 0);
  CHECK(empty_map.data() == nullptr);

  CHECK_THROWS_AS(MappedFile::map(tmp.file("absent.bin").string()), IoError);
}
