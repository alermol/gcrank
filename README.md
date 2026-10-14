# gcrank

Standalone C++17 library and CLI that answer GC-content queries over a FASTA
file (plain or gzip) in **O(1)** time, regardless of how large the window or
the genome is: a 250 Mb chromosome costs the same as a single base — two rank
lookups and a subtraction. Queries are answered from rank indices over two
bitvectors (GC and valid bases), so `N` and every IUPAC ambiguity code are
handled without special cases.

* **`gc_fraction()` is identical to `bedtools nuc`'s `pct_gc`** — same counts,
  same single-precision arithmetic, same printed digits; see
  [Comparison with `bedtools nuc`](#comparison-with-bedtools-nuc).
* Library `libgcrank`: C++17, standard library + zlib only (CRC-32, SHA-256
  and JSON are implemented in-tree); `main()` exists only in the CLI.
* CLI `gcrank-build`: builds an index from a plain or gzip FASTA (format
  sniffed by magic bytes, never by extension).
* Parallel builds that are **byte-identical for any thread count and chunk
  size**.
* `mmap` loading, 64-bit offsets throughout, versioned and checksummed on-disk
  format.
* **Tested on Linux** (GCC): Release, ASan+UBSan and ThreadSanitizer runs in
  CI. The Windows and macOS code paths are untested.

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build                                 # unit + CLI tests

./build/gcrank-build --fasta genome.fa.gz --out indexes --verbose
```

An index is written atomically to `<out>/<name>.gcrank/` (temp directory →
checksum verification → rename) and is never overwritten without
`--overwrite`.

Options: `--fasta <path>`, `--out <dir>`, `--name <name>` (default: input file
name minus the compression suffix and `.fa`/`.fasta`/`.fna`), `--chunk-size <N>`
(bases per parallel work unit, default 100,000,000), `--threads <N>` (default
`min(8, CPUs)`), `--overwrite`, `--verbose`, `--version`, `--help`.
Exit codes: 0 success, 1 usage, 2 input, 3 output, 4 internal.

## Library

```cpp
#include <cstdio>
#include "gcrank/gc_rank.h"

int main() {
  gcrank::GcRank index;

  gcrank::BuildOptions build;
  build.threads = 8;
  index.build("genome.fa.gz", build);         // parallel build
  index.save("indexes", "genome");            // -> indexes/genome.gcrank/

  gcrank::OpenOptions open;
  open.mmap = true;                           // map instead of reading
  gcrank::GcRank loaded;
  loaded.load("indexes/genome.gcrank", open); // falls back to reading on failure

  const uint64_t total = loaded.total_bases();
  std::printf("genome GC:  %.6f\n",
              loaded.gc_fraction(0, total, remove_ambiguous=false));
  std::printf("[1M, 2M):   %.6f\n",
              loaded.gc_fraction(1000000, 2000000, remove_ambiguous=true));
}
```

Consume the installed package with CMake:

```cmake
find_package(gcrank 1.0 REQUIRED)
target_link_libraries(app PRIVATE gcrank::gcrank)
```

### Query semantics

Windows are half-open `[start, end)`, 0-based, over the concatenated contigs
(per-contig offsets live in `contigs.tsv`).

| Call | Meaning |
| --- | --- |
| `gc_count(start, end)` | unambiguous `G`/`C` bases |
| `valid_count(start, end)` | unambiguous `A`/`C`/`G`/`T` bases |
| `gc_fraction(s, e, false)` | `GC / (e - s)` — the `bedtools nuc` convention |
| `gc_fraction(s, e, true)` | `GC / valid_count` — ambiguous bases removed |

`gc_fraction` converts the integer counts to `float` and divides in single
precision exactly like `bedtools nuc`, so `"%.6f"` output is bit-identical to
its `pct_gc` column; `gc_count` and `valid_count` stay exact integers.
`gc_fraction(..., true)` throws `QueryError` when the window holds no `A/C/G/T`
base, and zero-length windows are rejected either way. All errors derive from
`gcrank::Error` (`InputError`, `OutputError`, `FormatError`, `QueryError`,
`UsageError`, `IoError`).

### Base vocabulary

Accepted in either case: `A C G T U`, `N`, and the IUPAC codes
`R Y S W K M B D H V`. Every accepted byte counts towards `total_bases`, contig
lengths and `contigs.tsv`; the two bitvector bits decide what is counted where:

| Input | valid | GC | Notes |
| --- | --- | --- | --- |
| `A C G T` | yes | `G`/`C` only | unambiguous DNA |
| `U` | yes | no | folded onto `T`, counted the same way as `bedtools nuc` |
| `N` and IUPAC codes | no | no | possible `G`/`C` is never counted, so `gc_count` is a lower bound |

Anything else — including alignment gaps `-`, `.`, `*` — is rejected with
`InputError` before any index is written. The authoritative table is in
`include/gcrank/format.h`.

## Comparison with `bedtools nuc`

**GC-content values are identical.** `gc_fraction(start, end)` (the default,
`remove_ambiguous = false`) reproduces `bedtools nuc`'s `pct_gc` exactly, not
approximately: both tools count only unambiguous `G`/`C`, keep `N` and every
ambiguity code in the denominator (`seq_len`), fold `U` onto `T`, and convert
the counts to `float` before dividing — so the values match at all six printed
decimals with no tolerance.

| `bedtools nuc` | gcrank |
| --- | --- |
| `pct_gc` | `gc_fraction(start, end)` |
| `num_G + num_C` | `gc_count(start, end)` |
| `num_A + num_C + num_G + num_T` | `valid_count(start, end)` (`U` counted as `T`) |
| `seq_len` | `end - start` |

Verified against bedtools v2.31.1:

* `N` stays in the denominator: `AACCGGTTNNNN` → `0.333333` (4 of 12, not 4 of 8).
* Ambiguity codes are never GC and stay in the denominator: `AACCGGTTRRRR` →
  `0.333333`; `SSSSSSSSSS` → `0.000000` (`S` is "G or C", yet counts as other).
* Lowercase counts as a base: `aaccggtt` → `0.500000`.
* `U` is counted as `T`: `AACCGGTTUUUU` → `0.333333`.
* Rounding is shared: 20,000 natural windows (500–800 kb, with `N`) →
  20,000/20,000 identical; 16 adversarial windows with counts above 2^24 →
  16/16 identical (a double-precision quotient would differ on all 16). The
  counts agree exactly too (e.g. 273,711 GC over 603,818 bases prints
  `0.453301` in both tools; double rounding would give `0.453300`), and on a
  14.5 Gbp assembly the 21,000-window sums equalled the index totals
  (`GC = 6,695,634,623`, `ACGT = 14,509,208,007`).

Deliberate differences, none of them in `pct_gc`:

* `remove_ambiguous = true` (`GC / (A+C+G+T)`) has no `bedtools nuc` column.
* Alignment gaps (`-`, `.`, `*`) are rejected with `InputError` instead of
  being silently counted as `num_oth`.
* Zero-length windows and `start > end` throw `QueryError`; `bedtools nuc`
  skips the former (with a "Skipping." note) and errors on the latter.

The integer columns are exact, so prefer them when writing checks:

```cpp
// For each row of `bedtools nuc -fi genome.fa -bed windows.bed`
// (chrom, start, end, pct_at, pct_gc, num_A, num_C, num_G, num_T, num_N, num_oth, seq_len):
const gcrank::ContigInfo* contig = index.contig(chrom);
const uint64_t begin = contig->offset + start;   // index offsets are global
const uint64_t finish = contig->offset + end;
bool ok = (num_g + num_c == index.gc_count(begin, finish)) &&
          (num_a + num_c + num_g + num_t == index.valid_count(begin, finish)) &&
          (seq_len == end - start);
```

## Index format

```
<out>/<NAME>.gcrank/
├── manifest.json        canonical JSON: versions, source, checksums, semantics
├── contigs.tsv          <name>\t<length>\t<offset>, LF, no header
├── genome.bv            GC bitvector (1 = unambiguous G/C)
├── genome.rb            rank index for genome.bv
├── valid.bv             valid-base bitvector (1 = unambiguous A/C/G/T/U)
├── valid.rb             rank index for valid.bv
└── *.crc32              sidecar checksums for every data file
```

* Bit `i` is bit `(7 - i%8)` of byte `i/8` (MSB0); bytes are grouped into
  64-bit words with the first byte least significant — fixed by the format and
  host-independent.
* `.bv` files are zero-padded to whole 4096-bit superblocks. `.rb` files hold
  `ceil(n/4096) + 1` absolute ranks (`uint64`) followed by `ceil(n/512) + 1`
  ranks within their superblock (`uint16`), in the manifest's `endianness`;
  readers byte-swap when the host differs, so indexes move between hosts. Rank
  support costs 3 bytes per 512-bit block, 4.7% of a bitvector.
* `manifest.json` is the source of truth for the logical size and carries a
  SHA-256 of the source plus size/CRC-32 for every file; it advertises
  `remove_ambiguous_supported`. The format is documented in
  `include/gcrank/format.h`.

## Reproducibility, memory, performance

* Same input, options and `SOURCE_DATE_EPOCH` → byte-identical index,
  independent of the thread count and `--chunk-size` (both are scheduling
  knobs only). The manifest's `created` field is the only wall-clock read in
  the build path; enforced by `test_parallel.cpp`.
* Build peak memory ≈ `1.25 * chunk_bits * (threads + 2)` bytes, plus
  `2 * N/8` for the finished bitvectors and 4.7% of that for the rank indices.
  Lower `--chunk-size` to reduce the footprint. With `mmap = true`, resident
  memory is the rank arrays plus whatever the OS pages in; loading a 25 MB
  index takes ~0.15 s.
* Queries are rank arithmetic: ~230 ns per random window on a 100 Mb genome,
  and a whole-genome query costs ~290 ns — window length is irrelevant.

## Testing

```sh
ctest --test-dir build --output-on-failure      # doctest unit tests + CLI test
cmake -S . -B build-tsan -DGCRANK_TSAN=ON       # ThreadSanitizer build
```

The suite covers the rank index at every block boundary, FASTA parsing (plain,
CRLF, gzip, long lines, invalid bytes, every IUPAC code), `remove_ambiguous`
semantics, manifest round-trips, mmap vs heap equivalence, thread/chunk
invariance, and a committed golden index that is rebuilt and compared byte for
byte. Reference values come from `data-raw/make_expected_gc.py`, which never
calls gcrank; see `data-raw/README.md` for the regeneration policy.

gcrank is **tested on Linux** (GCC). CI
([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs the full suite in
three configurations on every push and pull request: Release (static and
shared), ASan + UBSan, and ThreadSanitizer. The Windows and macOS code paths
exist but are not covered by CI and must be treated as untested.

## Layout

```
include/gcrank/   public headers (API documentation in the headers)
src/              library implementation
cli/              gcrank-build
tests/            doctest suite, fixtures, golden index
data-raw/         fixture/reference generation scripts
```

MIT license — see `LICENSE`.
