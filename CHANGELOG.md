# Changelog

All notable changes to gcrank are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project uses
[Semantic Versioning](https://semver.org/). The on-disk format carries its own
`format_version`, which is versioned independently of the package.

## [0.2.0] - 2026-10-25

### Added

- `GcRank::stats()`, which returns the window length, GC count, valid count and
  GC fraction from one contig lookup and one rank pass per bitvector. Undefined
  fractions (empty window, or no A/C/G/T base with `remove_ambiguous`) are
  reported as NaN instead of throwing.

### Changed

- Faster index builds: SHA-256 uses the x86 SHA extensions when the CPU has
  them (runtime dispatched, scalar fallback), hashing overlaps parsing on a
  helper thread, and bases are packed eight at a time with a lookup table.
- Faster queries: the CLI resolves each window once through `stats()`, formats
  rows with `std::to_chars`, and uses a larger stdout buffer.

### Removed

- Unused internal helpers: `detail::println()`, `CRC32::reset()`,
  `ThreadPool::size()`, `PlainInputStream(FILE*, std::string)`, and the
  `close_enough()` test helper.
- Includes that had no direct users.

No on-disk format changes.

## [0.1.0] - 2026-10-16

First release.
