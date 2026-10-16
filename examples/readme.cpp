// The library example from README.md, verbatim apart from this comment and the
// blank-line/indent layout.
//
// The test build compiles this target, so the documented example is guaranteed
// to compile against the public API. It is never run: "genome.fa.gz" does not
// exist in the source tree.
#include <cstdio>

#include "gcrank/gc_rank.h"

int main() {
  gcrank::GcRank index;

  gcrank::BuildOptions build;
  build.threads = 8;
  index.build("genome.fa.gz", build);          // parallel build
  index.save("indexes", "genome");             // -> indexes/genome.gcrank/

  gcrank::OpenOptions open;
  open.mmap = true;                            // map instead of reading
  gcrank::GcRank loaded;
  loaded.load("indexes/genome.gcrank", open);  // falls back to reading on failure

  const uint64_t total = loaded.total_bases();
  std::printf("genome GC:  %.6f\n",
              loaded.gc_fraction(0, total, /*remove_ambiguous=*/false));
  std::printf("[1M, 2M):   %.6f\n",
              loaded.gc_fraction(1000000, 2000000, /*remove_ambiguous=*/true));
}
