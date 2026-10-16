# Changelog

All notable changes to gcrank are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/). The on-disk format carries its own
`format_version`, which is versioned independently of the package.

## [0.1.0] - 2026-10-16

First release.

### Added

- `libgcrank`: a C++ library (C++20 to consume the installed headers) with
  `GcRank::build`/`save`/`load`, parallel FASTA ingestion, and O(1) window
  queries through `gc_count`, `valid_count` and `gc_fraction`, in global and
  contig-relative coordinates. Only the standard library and zlib are needed:
  CRC-32 uses zlib, SHA-256 and canonical JSON are implemented in-tree.
- `gcrank-build`: builds an index from a plain or gzip FASTA file (the format
  is sniffed from the magic bytes, never the extension; `--fasta -` reads
  stdin and pipes), with atomic writes, `--overwrite`, `--threads`,
  `--chunk-size` and `--verbose`.
- `gcrank-query`: answers BED intervals or genome-wide tilings from an index
  (`--bed`, `--window`/`--step`, `--remove-ambiguous`, `--mmap`,
  `--no-validate-crc`), with `bedtools`-identical `pct_gc` output.
- On-disk format `0.1.0`: `manifest.json` (canonical JSON with the source
  SHA-256, sizes and CRC-32 for every file, and the `total_gc` / `total_valid`
  population counts), `contigs.tsv`, two bitvectors (`genome.bv`: unambiguous
  G/C; `valid.bv`: unambiguous A/C/G/T) with two-level rank indices, and
  sidecar checksums.
- Queries are half-open, 0-based windows with two denominators: GC / window
  length (default) and GC / (A+C+G+T) with `remove_ambiguous = true`.
- Base vocabulary: `A C G T U`, `N` and the IUPAC ambiguity codes
  `R Y S W K M B D H V`, in either case. Only unambiguous G/C counts as GC;
  `U` is folded onto `T`. Alignment gaps (`-`, `.`, `*`) are rejected with
  `InputError` instead of being silently miscounted, and duplicate or empty
  contig names are rejected as well.
- `gc_fraction` is bit-identical to `bedtools nuc`'s `pct_gc`: the integer
  counts are converted to `float` and divided in single precision. The parity
  is re-verified against the real `bedtools nuc` on an exhaustive fixture
  window set and on adversarial counts above 2^24
  (`data-raw/check_bedtools_parity.py`, run by CTest as the `bedtools` test).
- Parallel builds that are byte-identical for any thread count and
  `--chunk-size`; the source is hashed while it is parsed (one pass, no
  hash/parse mismatch window), and `SOURCE_DATE_EPOCH` pins the manifest
  timestamp, making builds reproducible for the same input path.
- `mmap` loading with a transparent fallback to reading. The mapping is
  advisory (`MADV_RANDOM`), and the SIGBUS-on-truncation caveat is documented.
- An optional log sink on `BuildOptions`/`OpenOptions`: the library reports
  progress and warnings through it and never writes to stderr itself; the CLIs
  install the stderr sink.
- CMake package export (`find_package(gcrank)`, `gcrank::gcrank`); only the
  four public headers are installed, everything else lives under
  `include/gcrank/detail/`.
- Test suite (doctest) covering rank-index boundaries and validation, FASTA
  parsing (plain, CRLF, gzip, truncated gzip, multi-member gzip, unwrapped
  single-line records, invalid bytes, every IUPAC code), query semantics,
  contig-relative queries, manifest round-trips, mmap vs heap equivalence,
  thread/chunk invariance, concurrent saves, a committed golden index rebuilt
  byte for byte, the documented bit layout (`data-raw/check_spec.py`), the
  bedtools parity check, and an install + `find_package` consumer test. CI runs
  the suite on Linux in Release (`-Werror`, static and shared), Clang 18 with
  libstdc++ 14, ASan + UBSan and ThreadSanitizer configurations.

### Notes

- A truncation at exactly a member boundary of a multi-member gzip file cannot
  be detected, because gzip records no member count.
- A `save` that crashes between its two renames can leave a
  `<name>.gcrank.tmp.*` or `<name>.gcrank.old.*` directory behind. Unique
  names guarantee that a later save never reuses or deletes one; removing such
  leftovers is a manual step.
