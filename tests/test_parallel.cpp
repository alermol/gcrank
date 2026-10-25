#include <doctest.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "gcrank/error.h"
#include "gcrank/gc_rank.h"
#include "gcrank/detail/thread_pool.h"
#include "test_util.h"

using namespace gcrank;

namespace {

namespace fs = std::filesystem;

// Same source and same SOURCE_DATE_EPOCH: only the thread count differs, which
// must not change a single byte of the output.
std::string build_and_diff(gctest::TempDir& tmp, unsigned threads, uint64_t chunk_bases) {
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", "1700000000");
  GcRank index;
  BuildOptions options;
  options.threads = threads;
  options.chunk_bits = chunk_bases;
  index.build(gctest::data_file("test_small.fa").string(), options);
  const std::string name = "t" + std::to_string(threads) + "c" + std::to_string(chunk_bases);
  index.save(tmp.path().string(), name);
  return (tmp.file(name + ".gcrank")).string();
}

}  // namespace

TEST_CASE("parallel: identical bytes for 1, 2, 4, 8 and 16 threads") {
  gctest::TempDir tmp("parallel-threads");
  std::vector<std::string> dirs;
  for (unsigned threads : {1u, 2u, 3u, 4u, 8u, 16u}) {
    CAPTURE(threads);
    dirs.push_back(build_and_diff(tmp, threads, 997));
  }
  for (std::size_t i = 1; i < dirs.size(); ++i) {
    CAPTURE(i);
    CHECK(gctest::diff_dirs(dirs[0], dirs[i]).empty());
  }
}

TEST_CASE("parallel: the chunk size does not change the output") {
  gctest::TempDir tmp("parallel-chunks");
  std::vector<std::string> dirs;
  for (uint64_t chunk : {64ull, 997ull, 4096ull, 100000000ull}) {
    CAPTURE(chunk);
    dirs.push_back(build_and_diff(tmp, 4, chunk));
  }
  for (std::size_t i = 1; i < dirs.size(); ++i) {
    CAPTURE(i);
    CHECK(gctest::diff_dirs(dirs[0], dirs[i]).empty());
  }
}

TEST_CASE("parallel: a larger synthetic genome is consistent across thread counts") {
  gctest::TempDir tmp("parallel-large");
  const fs::path fasta = tmp.file("synthetic.fa");
  const uint64_t contig_bases = 250000;
  const uint64_t total_bases = contig_bases * 9;
  gctest::write_synthetic_fasta(fasta, total_bases, contig_bases, /*n_every=*/13);

  std::string first_dir;
  for (unsigned threads : {1u, 2u, 4u, 8u}) {
    CAPTURE(threads);
    const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", "1700000000");
    GcRank index;
    BuildOptions options;
    options.threads = threads;
    options.chunk_bits = 50000;  // ~45 blocks
    index.build(fasta.string(), options);
    CHECK(index.total_bases() == total_bases);
    CHECK(index.contigs().size() == 9);

    // Spot-check against an independent count of the generated file.
    const gctest::NaiveFasta reference = gctest::parse_fasta_plain(fasta);
    CHECK(index.gc_count(0, total_bases) == gctest::naive_gc_count(reference, 0, total_bases));
    CHECK(index.valid_count(0, total_bases) ==
          gctest::naive_valid_count(reference, 0, total_bases));
    CHECK(index.gc_count(12345, 987654) == gctest::naive_gc_count(reference, 12345, 987654));
    CHECK(index.valid_count(12345, 987654) ==
          gctest::naive_valid_count(reference, 12345, 987654));

    const std::string name = "large" + std::to_string(threads);
    index.save(tmp.path().string(), name);
    const std::string dir = tmp.file(name + ".gcrank").string();
    if (first_dir.empty()) {
      first_dir = dir;
    } else {
      CHECK(gctest::diff_dirs(first_dir, dir).empty());
    }
  }
}

TEST_CASE("parallel: task exceptions keep their type through wait() and stop()") {
  // A failing task must surface as its original exception type on the thread
  // that owns the pool, not as the std::exception_ptr wrapper.
  ThreadPool pool(1);
  pool.submit([] { throw InputError("worker blew up"); });
  CHECK_THROWS_AS(pool.wait(), InputError);

  ThreadPool second(2);
  second.submit([] { throw QueryError("bad window"); });
  CHECK_THROWS_AS(second.stop(), QueryError);
}

TEST_CASE("parallel: threads = 0 means auto and still matches the single-thread output") {
  gctest::TempDir tmp("parallel-zero");
  const gctest::EnvGuard epoch("SOURCE_DATE_EPOCH", "1700000000");
  GcRank index;
  BuildOptions options;
  options.threads = 0;
  options.chunk_bits = 1000;
  CHECK_NOTHROW(index.build(gctest::data_file("test_small.fa").string(), options));
  index.save(tmp.path().string(), "zero");
  const std::string one = build_and_diff(tmp, 1, 1000);
  CHECK(gctest::diff_dirs(tmp.file("zero.gcrank").string(), one).empty());
}
