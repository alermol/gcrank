// Consumer smoke test for the installed package: builds with C++20 against
// find_package(gcrank) and exercises build/save/load and the query API.
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "gcrank/gc_rank.h"

int main() {
  const std::filesystem::path dir = std::filesystem::temp_directory_path() / "gcrank-consumer";
  std::error_code ec;
  std::filesystem::remove_all(dir, ec);
  std::filesystem::create_directories(dir, ec);
  if (ec) {
    std::fprintf(stderr, "cannot create %s\n", dir.string().c_str());
    return 1;
  }

  const std::filesystem::path fasta = dir / "tiny.fa";
  {
    std::ofstream out(fasta);
    out << ">chr1\nACGTACGTNN\n>chr2\nGGCC\n";
  }

  gcrank::GcRank index;
  index.build(fasta.string());
  index.save(dir.string(), "tiny");

  gcrank::GcRank loaded;
  loaded.load((dir / "tiny.gcrank").string());
  if (loaded.total_bases() != 14) return 2;
  if (loaded.gc_count("chr1", 0, 10) != 4) return 3;
  if (loaded.valid_count("chr2", 0, 4) != 4) return 4;
  if (loaded.contig("chr1") == nullptr || loaded.contig("chr1")->offset != 0) return 5;

  std::filesystem::remove_all(dir, ec);
  std::printf("consumer ok: %llu bases, %.6f GC\n",
              static_cast<unsigned long long>(loaded.total_bases()),
              loaded.gc_fraction(0, loaded.total_bases(), false));
  return 0;
}
