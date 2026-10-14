# Changelog

All notable changes to gcrank are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html). The on-disk format
carries its own `format_version`, which is versioned independently of the
package.

## [Unreleased]

### Changed

- The library and CLI now require C++23 (GCC 14+, Clang 18+ with libstdc++ 14+,
  or MSVC 19.40+). The on-disk format, CLI output and query results are
  unchanged; the golden-index byte-for-byte tests still pass.
- Modernized the internals with the C++20/23 standard library: `std::print`
  and `std::format` for diagnostics and canonical JSON, `std::span` for
  bitvector rank/word views and packed blocks, `<bit>` (`std::byteswap`,
  `std::endian`, `std::popcount`, `std::rotr`) replacing hand-rolled
  byte-swapping, the population count and the CRC-32 table (now computed at
  compile time), `std::move_only_function` for thread-pool tasks,
  `std::jthread` with a stop token for the producer, `std::from_chars` for
  parsing and `std::format` for rendering, and range algorithms throughout.
  `BaseClass` is now a scoped enum and `BitVector::append_bits` takes a
  `std::span`.

## [1.0.0] - 2026-10-14

First release.

### Added

- `libgcrank`: a C++17 library with no dependencies beyond zlib and the
  standard library (CRC-32, SHA-256 and canonical JSON are implemented
  in-tree). `GcRank::build`/`save`/`load`, parallel FASTA ingestion, and O(1)
  window queries through `gc_count`, `valid_count` and `gc_fraction`.
- `gcrank-build`: a CLI that builds an index from a plain or gzip FASTA file
  (the format is sniffed from the magic bytes, never the extension), with
  atomic writes, `--overwrite`, `--threads`, `--chunk-size` and `--verbose`.
- On-disk format `1.0.0`: `manifest.json` (canonical JSON with the source
  SHA-256, sizes and CRC-32 for every file), `contigs.tsv`, two bitvectors
  (`genome.bv`: unambiguous G/C; `valid.bv`: A/C/G/T) with two-level rank
  indices, and sidecar checksums.
- O(1) queries over global, 0-based half-open windows, with two denominators:
  GC / window length (default) and GC / (A+C+G+T) with
  `remove_ambiguous = true`.
- Base vocabulary: `A C G T U`, `N` and the IUPAC ambiguity codes
  `R Y S W K M B D H V`, in either case. Only unambiguous G/C counts as GC;
  `U` is folded onto `T`. Alignment gaps (`-`, `.`, `*`) are rejected with
  `InputError` instead of being silently miscounted.
- `gc_fraction` is bit-identical to `bedtools nuc`'s `pct_gc`: the integer
  counts are converted to `float` and divided in single precision, verified on
  20,000 natural windows and 16 adversarial windows with counts above 2^24.
- Parallel build whose output is byte-identical for any thread count and
  `--chunk-size`; `SOURCE_DATE_EPOCH` pins the manifest timestamp, making
  builds fully reproducible.
- `mmap` loading for large indexes, with a transparent fallback to reading.
- CMake package export (`find_package(gcrank)`, `gcrank::gcrank`) for
  downstream consumers.
- Test suite (doctest: 58 cases, 92,793 assertions) covering rank-index
  boundaries, FASTA parsing (plain, CRLF, gzip, long lines, invalid bytes,
  every IUPAC code), query semantics, manifest round-trips, mmap vs heap
  equivalence, thread/chunk invariance, and a committed golden index rebuilt
  and compared byte for byte. CI runs the suite on Linux with Release
  (static and shared), ASan + UBSan and ThreadSanitizer builds.
