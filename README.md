# gcrank

Standalone C++23 library and CLI that answer GC-content queries over a FASTA
file (plain or gzip) in **O(1)** time, regardless of how large the window or
the genome is: a 250 Mb chromosome costs the same as a single base — two rank
lookups and a subtraction. Queries are answered from rank indices over two
bitvectors (GC and valid bases), so `N` and every IUPAC ambiguity code are
handled without special cases.

* **`gc_fraction()` is identical to `bedtools nuc`'s `pct_gc`** — same counts,
  same single-precision arithmetic, same printed digits; see
  [Comparison with `bedtools nuc`](#comparison-with-bedtools-nuc).
* Library `libgcrank`: C++23 (C++20 to consume the installed headers),
  standard library + zlib only (CRC-32 via zlib, SHA-256 and JSON implemented
  in-tree); `main()` exists only in the CLIs.
* CLI `gcrank-build`: builds an index from a plain or gzip FASTA (format
  sniffed by magic bytes, never by extension), including pipes and
  `/dev/stdin` via `--fasta -`.
* CLI `gcrank-query`: answers BED intervals or whole-contig tilings from an
  index, with `bedtools`-identical `pct_gc` output; no C++ needed to query.
* Contig-relative queries (`gc_count("chr1", start, end)`), so BED coordinates
  can be passed through unchanged instead of doing offset arithmetic.
* Parallel builds that are **byte-identical for any thread count and chunk
  size**, with the source hashed while it is parsed (one pass, no
  hash/parse mismatch window).
* `mmap` loading, 64-bit offsets throughout, versioned and checksummed on-disk
  format (sidecar checksums and rank totals are cross-checked on every load).
* **Tested on Linux** (GCC 14 and Clang 18): Release, ASan+UBSan and
  ThreadSanitizer runs in CI. The Windows and macOS code paths are untested.

## Requirements

* To build: a C++23 compiler (CI uses GCC 14 and Clang 18 with libstdc++ 14;
  GCC 13 is expected to work as well, untested). The build uses `std::format`,
  `std::span`, `<bit>` (`std::byteswap`, `std::endian`, `std::popcount`,
  `std::rotr`), `std::move_only_function`, `std::jthread`, `std::from_chars`
  and the C++20/23 range algorithms.
* To consume the installed library: C++20 (`std::span`, `std::string_view`);
  only the four public headers (`error.h`, `format.h`, `gc_rank.h`,
  `manifest.h`) are installed, everything else lives under `detail/`.
* CMake 3.20+ and zlib (the only third-party dependency).

## Quick start

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build                                 # unit + CLI + spec + consumer

./build/gcrank-build --fasta genome.fa.gz --out indexes --verbose
echo -e 'chr1\t1000\t2000' | ./build/gcrank-query --index indexes/genome.gcrank --bed -
```

An index is written atomically to `<out>/<name>.gcrank/` (uniquely named temp
directory → checksum verification → `fsync` → rename) and is never overwritten
without `--overwrite`. With `--overwrite`, the new index is fully written and
verified before the existing one is moved aside, and it is put back if the
swap fails.

Options: `--fasta <path>` (or `-` for stdin), `--out <dir>`, `--name <name>`
(default: the input file name with a compression or FASTA suffix removed, or
`stdin` for `-`), `--chunk-size <N>` (bases per parallel work unit, default
16,000,000, at most 1 GiB), `--threads <N>` (default `min(8, CPUs)`, at most
1024), `--overwrite`, `--verbose`, `--version`, `--help`.
Exit codes: 0 success, 1 usage, 2 input, 3 output, 4 internal.
An obviously oversized `--chunk-size × --threads` combination warns before
allocating.

## Querying an index: `gcrank-query`

```sh
# BED intervals (chrom, start, end; 0-based half-open, like bedtools)
gcrank-query --index indexes/genome.gcrank --bed windows.bed > windows.tsv

# the same intervals with ambiguous bases removed from the denominator
gcrank-query --index indexes/genome.gcrank --bed windows.bed --remove-ambiguous

# genome-wide tiling, no BED file needed
gcrank-query --index indexes/genome.gcrank --window 10000 --step 5000 > tiles.tsv
```

Output is tab separated with a `#` header, one row per window:

```
#chrom  start  end  pct_gc  num_gc  num_valid  seq_len
```

`pct_gc` is printed with `%.6f` and is bit-identical to `bedtools nuc`'s
`pct_gc` for the same interval (`num_gc` is `num_C + num_G`; `num_valid` is
`num_A + num_C + num_G + num_T`, with `U` counted as `T`). `--remove-ambiguous`
changes the denominator to `num_valid` and prints `nan` where there is no
`A/C/G/T` base. A window that leaves its contig, or names an unknown contig, is
reported as an input error with the BED line number instead of silently
crossing into the next contig. `--bed -` reads BED from stdin; `--mmap` and
`--no-validate-crc` trade start-up work for speed on huge indexes.

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
              loaded.gc_fraction(0, total, /*remove_ambiguous=*/false));
  std::printf("[1M, 2M):   %.6f\n",
              loaded.gc_fraction(1000000, 2000000, /*remove_ambiguous=*/true));
}
```

Consume the installed package with CMake:

```cmake
find_package(gcrank 0.3 REQUIRED)
target_link_libraries(app PRIVATE gcrank::gcrank)
```

### Query semantics

Windows are half-open `[start, end)`, 0-based. Every query exists twice: with
global coordinates over the concatenated contigs (per-contig offsets live in
`contigs.tsv`), and with contig-relative coordinates
(`index.gc_count("chr1", start, end)`), which reject windows that leave the
named contig with `QueryError`. Prefer the contig-relative form: BED files are
contig-relative, and using global coordinates without adding the offset
returns the neighbouring contig's answer.

| Call | Meaning |
| --- | --- |
| `gc_count(start, end)` | unambiguous `G`/`C` bases |
| `valid_count(start, end)` | unambiguous `A`/`C`/`G`/`T` bases |
| `gc_fraction(s, e, false)` | `GC / (e - s)` — the `bedtools nuc` convention |
| `gc_fraction(s, e, true)` | `GC / valid_count` — ambiguous bases removed |
| `gc_count("chr1", s, e)` | the same, relative to the contig start |
| `stats("chr1", s, e, rm)` | `{length, gc, valid, gc_fraction}` in one pass |

`gc_fraction` converts the integer counts to `float` and divides in single
precision exactly like `bedtools nuc`, so `"%.6f"` output is bit-identical to
its `pct_gc` column; `gc_count` and `valid_count` stay exact integers.
`gc_fraction(..., true)` throws `QueryError` when the window holds no `A/C/G/T`
base, and zero-length windows are rejected either way. `stats()` is the
efficient choice when several columns are needed: it resolves the contig once,
reads each rank index once, and reports an undefined fraction as NaN instead
of throwing. All errors derive from
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
bool ok = (num_g + num_c == index.gc_count(chrom, start, end)) &&
          (num_a + num_c + num_g + num_t == index.valid_count(chrom, start, end)) &&
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

* Bit `i` is bit `63 - i%64` of 64-bit little-endian word `i/64`, i.e. bit
  `(7 - i%8)` of byte `8*(i/64) + 7 - (i%64)/8` of the file (MSB0 within the
  word) — fixed by the format and host-independent. `data-raw/check_spec.py`
  (run by CTest as the `spec` test) decodes the golden index with an
  independent reader that implements exactly this rule.
* `.bv` files are zero-padded to whole 4096-bit superblocks. `.rb` files hold
  `ceil(n/4096) + 1` absolute ranks (`uint64`) followed by `ceil(n/512) + 1`
  ranks within their superblock (`uint16`). They are always little-endian; the
  manifest's `endianness` records the byte order, and the reader swaps when an
  index is stored in the opposite order. Rank support costs 3 bytes per
  512-bit block, 4.7% of a bitvector.
* `manifest.json` is the source of truth for the logical size; it carries a
  SHA-256 of the parsed source stream (identical for a FASTA and its gzip
  form, since the hash is computed while the input is read), size/CRC-32 for
  every file, the `index.total_gc` / `index.total_valid` population counts, and
  the advertised `remove_ambiguous_supported`. Loading always cross-checks the
  sidecar `.crc32` files and the rank totals against the manifest; the full
  data-file CRC pass runs unless `validate_crc = false`. The format is
  documented in `include/gcrank/format.h`.

## Reproducibility, memory, performance

* Same input, options, `SOURCE_DATE_EPOCH` **and path string** → byte-identical
  index, independent of the thread count and `--chunk-size` (both are
  scheduling knobs only). The manifest records the source path exactly as
  given, so a relative and an absolute build of the same file differ in that
  field; use the same path when comparing. The manifest's `created` field is
  the only wall-clock read in the build path; enforced by `test_parallel.cpp`.
* Build peak memory ≈ `1.25 * chunk_bits * (threads + 2)` bytes, plus
  `2 * N/8` for the finished bitvectors, 4.7% of that for the rank indices,
  and about 1 MB of hashing buffers. The default chunk is 16 M bases, so a
  300 Mb genome builds in about 190 MB of RSS regardless of thread count — and
  regardless of whether the FASTA is wrapped at 60 columns or one single line
  per contig, because sequence lines are never materialised.
* Hashing is not on the critical path: the SHA-256 helper thread consumes the
  parser's read buffers while the parallel workers pack blocks. SHA-NI is used
  when available (runtime detected, scalar fallback elsewhere). On a 160 Mbp
  60-column FASTA with an AMD EPYC 7742, a Release build takes ~0.8 s end to
  end for 1–8 threads; the remaining limit is the serial FASTA parse.
* Queries are rank arithmetic and independent of the window length: on that
  160 Mbp index, a random window costs ~0.3 µs per count (~0.65 µs for a full
  `stats()` row: `gc`, `valid`, `length` and the fraction) and a whole-genome
  query ~50 ns. Loading with `mmap = true` (`validate_crc = false`) maps the
  index instead of reading it; a 160,000-window tiling query takes ~0.19 s end
  to end, and 200,000 BED windows take ~0.32 s, including formatting.

## Testing

```sh
ctest --test-dir build --output-on-failure   # unit + CLI + spec + bedtools + consumer
cmake -S . -B build-tsan -DGCRANK_TSAN=ON    # ThreadSanitizer build
```

CI installs `bedtools` for the Release and Clang jobs, so the `bedtools` test
runs there; in configurations without it the test skips (CTest reports the
skip) and the claim can be re-checked by hand with
`python3 data-raw/check_bedtools_parity.py`.

The suite covers the rank index at every block boundary, FASTA parsing (plain,
CRLF, gzip, truncated gzip, multi-member gzip, long and single-line records,
invalid bytes, every IUPAC code), `remove_ambiguous` semantics, manifest
round-trips, mmap vs heap equivalence, thread/chunk invariance, the contig-
relative API, and a committed golden index that is rebuilt and compared byte
for byte. Reference values come from `data-raw/make_expected_gc.py`, which
never calls gcrank, and the documented bit layout is checked by
`data-raw/check_spec.py`; see `data-raw/README.md` for the regeneration policy.
The `consumer` test installs into a scratch prefix and builds a C++20 project
against `find_package(gcrank)`, so the exported target and the installed
headers are tested as a downstream user sees them.

gcrank is **tested on Linux** (GCC 14 and Clang 18). CI
([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs the full suite on
every push and pull request: Release with `-Werror` (static and shared),
Clang 18 with libstdc++ 14, ASan + UBSan, and ThreadSanitizer. The Windows and
macOS code paths exist but are not covered by CI and must be treated as
untested (AppleClang's libc++ does not implement `std::format` yet).

## Layout

```
include/gcrank/       public headers (error, format, gc_rank, manifest)
include/gcrank/detail internal headers (installed only in the source tree)
src/                  library implementation
cli/                  gcrank-build, gcrank-query
tests/                doctest suite, fixtures, golden index, consumer test
data-raw/             fixture/reference generation scripts and the spec checker
```

MIT license — see `LICENSE`.
