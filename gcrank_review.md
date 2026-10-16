# gcrank 2.0.0 — code review and improvement brief

**Audience:** an AI coding agent that will implement these changes.
**Scope reviewed:** every file in `gcrank-main/` except the vendored `tests/doctest.h`.
**How it was checked:** I read all the code, then built it (GCC 13 + a 10-line `<print>` shim, Debug with ASan+UBSan, plus Release). The existing suite passes. Every item marked **[verified]** was reproduced with a command or test program. The rest come from reading the code, and the reasoning is given.

## Ground rules for the implementing agent

1. **Do not break the on-disk format silently.** `tests/data/test_small.gcrank/` is a golden index that is rebuilt and compared byte-for-byte. If a change has to alter bytes (only P1-6 might), bump `kFormatVersion`, regenerate with `data-raw/make_golden_index.sh`, and record it in `CHANGELOG.md`.
2. **Add a failing test first** for every bug in P0/P1, then fix it. Most items below list the test they need.
3. Keep `gc_fraction()` bit-identical to `bedtools nuc` by default. This is the project's headline claim.
4. Run all three CI configurations (Release, ASan+UBSan, TSan) before you finish.

---

## P0 — Correctness bugs (fix first)

### P0-1. A truncated gzip input is silently accepted as a shorter genome  **[verified]**
- **Where:** `src/input_stream.cpp:60-72` (`GzipInputStream::read`).
- **Problem:** `gzread()` returns `0`, not `-1`, when the stream ends early. zlib deliberately leaves `Z_BUF_ERROR` ("unexpected end of file") out of its error return. The code treats `0` as a clean EOF.
- **Evidence:** I cut a 2 Mb `.fa.gz` in half. `gcrank-build` exited 0 and wrote an index with `total_bases: 999090` instead of 2,000,000. The same happened with the test fixture cut at 2500 bytes: exit 0, 5 of 8 contigs.
- **Impact:** A partial download or copy turns into a valid-looking, checksummed index of the wrong genome. This is the worst bug in the codebase.
- **Fix:** When `gzread` returns 0, call `gzerror(gz, &err)`. If `err != Z_OK`, throw `InputError` (this includes `Z_BUF_ERROR`). Also check the return value of `gzclose_r`. Optionally, at EOF, use `gzdirect()` and the gzip trailer to confirm the stream finished.
- **Test:** truncate `test_small.fa.gz` at several offsets and require `InputError` / CLI exit code 2. Also check that a valid multi-member gzip still parses.

### P0-2. `ThreadPool::wait()` / `stop()` throw the `std::exception_ptr` object instead of rethrowing it  **[verified]**
- **Where:** `src/thread_pool.cpp:70` and `:87` (`if (error) throw error;`).
- **Problem:** `throw error;` throws an object of type `std::exception_ptr`. No `catch (const std::exception&)` or `catch (const gcrank::Error&)` can catch it.
- **Evidence:** a task that throws `gcrank::InputError` reached the caller as a `std::exception_ptr`.
- **Impact:** `thread_pool.h` is an installed public header, so downstream users hit this. Inside gcrank, `BuildPipeline::shutdown()` swallows the exception with `catch (...)`, which hides the bug.
- **Fix:** `std::rethrow_exception(error);` in both places. Then see P2-1: the class probably should not be public at all.
- **Test:** submit a task that throws `InputError`, then `CHECK_THROWS_AS(pool.wait(), InputError)`.

### P0-3. `--overwrite` deletes the existing index before the new one is written (not atomic)
- **Where:** `src/gc_rank.cpp:569-575` (`GcRank::save`).
- **Problem:** `remove_all(final_dir)` runs **before** the temporary directory is written and verified. If anything fails afterwards (disk full, I/O error, a crash), the old index is gone and no new one exists. The README says "written atomically", which is false for overwrites.
- **Other problems in the same function:**
  - The temporary name is fixed (`<name>.gcrank.tmp`, line 563). Two concurrent saves of the same name delete each other's work (`remove_all(tmp_dir)` at 576-579).
  - There is no `fsync` of the files or the parent directory before `rename`. `verify_index_files` re-reads from the page cache, so it proves nothing about durability.
  - `fs::exists(final_dir, ec)` ignores `ec`.
- **Fix:**
  1. Write to a unique temp directory (for example `<name>.gcrank.tmp.<pid>.<random>`).
  2. Verify it, then `fsync` the files and the directory.
  3. If the target exists, rename it to `<name>.gcrank.old.<rand>`, rename tmp → final, then remove the old one. On Linux ≥ 3.15 you can use `renameat2(RENAME_EXCHANGE)` for a true swap.
  4. Remove only the temp directory this call created.
- **Test:** inject a failure after the old index exists (for example, a test hook or an output directory made read-only partway through). Assert the old index still loads.

### P0-4. Duplicate and empty contig names are accepted  **[verified]**
- **Where:** `src/fasta_reader.cpp:51-61` (`parse_header`). Nothing checks names later either.
- **Evidence:** the input `>chr1 / >chr1 / > / > ` built successfully. `contigs.tsv` ended up with two `chr1` rows and two rows with an empty name.
- **Impact:** `GcRank::contig(name)` silently returns the first match. Any BED-driven query against a duplicated name returns the wrong contig's GC.
- **Fix:** Reject an empty name with `InputError` and the line number. Reject duplicate names with `InputError`; checking at assembly time with an `unordered_set` is enough. `samtools faidx` also refuses duplicates.
- **Test:** add FASTA fixtures for both cases and expect `InputError` / exit code 2.

### P0-5. The format spec documents the wrong bit layout  **[verified]**
- **Where:** `include/gcrank/format.h:15-18`, `README.md:188-190`, `include/gcrank/bitvector.h:5-7`. These also contradict each other.
- **Evidence:** I decoded the golden `genome.bv` with a Python script and compared it to the FASTA:
  - The documented rule ("bit *i* = bit (7 − i%8) of byte i/8", i.e. plain MSB0 bytes) does **not** match.
  - The actual layout is: each 64-bit word holds bit *i* at position 63 − i%64, and the word is stored **little-endian**. So bit *i* is bit `(7 − i%8)` of byte `8·(i/64) + 7 − (i%64)/8`. The bytes inside each 8-byte group come out reversed relative to the doc.
  - `bitvector.h` says "the very first byte of the stream occupies the most significant byte of word 0", but `load_le64` puts it in the **least** significant byte.
- **Impact:** Anyone writing a reader in Python or R from the spec gets the wrong bits. The format is presented as a "public contract".
- **Fix:** Correct the documentation in all three places so it matches the bytes the code actually writes. Do not change the code, because the golden index and existing indexes depend on it. Add a doc-test: a tiny Python reader in `data-raw/` that decodes the golden index using only the documented rule and checks it against `expected_gc.tsv`.

### P0-6. Documentation says `U` is ambiguous; the code treats it as a valid base (T)
- **Where:** `format.h:22` ("Ambiguous input (N, U, R, …) sets neither bit"), `gc_rank.h:73-74` and `:77-78` (both list U among the excluded codes).
- **Truth:** `U` sets the valid bit (`format.h:174-175`, `default_semantics()`, the README table).
- **Fix:** Correct the three comments. Also fix the typo "multilbyte" at `format.h:25`.

### P0-7. The README library example does not compile
- **Where:** `README.md:74-76`. The calls `gc_fraction(0, total, remove_ambiguous=false)` use Python-style named arguments, which C++ does not have.
- **Fix:** use `/*remove_ambiguous=*/false`. Add a CI step that compiles the README snippet (extract it into `examples/readme.cpp` and build it).

### P0-8. CLI exits with code 4 after successfully writing an index for an empty FASTA with `--verbose`  **[verified]**
- **Where:** `cli/gcrank_build.cpp:192-194`. It calls `gc_fraction(0, total_bases())` with `total_bases() == 0`, which throws `QueryError`.
- **Evidence:** `gcrank-build --fasta empty.fa --out o --verbose` wrote `empty.gcrank` and then exited 4 with "error: gc_fraction: empty window".
- **Fix:** skip the fraction (or print `n/a`) when `total_bases() == 0`. Add a CLI test for it.

---

## P1 — Robustness, API design, performance

### P1-1. Queries use global coordinates only, and windows may cross contig boundaries  **[verified]**
- **Where:** `gc_rank.h` query API; `GcRank::check_range` (`gc_rank.cpp:634`).
- **Evidence:** `gc_fraction(50, 150)` over `chr_boundary_100` (length 100) followed by the next contig returned `0.540000` with no error. `bedtools nuc` would reject the equivalent BED interval.
- **Impact:** This is the most common real-world misuse: users take BED coordinates and forget to add `contig->offset`. The README even has users do that offset arithmetic by hand.
- **Fix:** Add contig-relative overloads, `gc_count(std::string_view contig, uint64_t start, uint64_t end)` and the same for `valid_count` and `gc_fraction`. They should check `end <= contig.length` and throw `QueryError` otherwise. Keep the global versions (perhaps rename them `*_global`) for advanced use.
- `contig(name)` is a linear scan (`gc_rank.cpp:682`). Build an `unordered_map<string, size_t>` at load and build time. Assemblies with 10⁵–10⁶ scaffolds and large BED files make the current version O(n·m).

### P1-2. No query tool, only a builder
- Users must write C++ to query an index. Add `gcrank-query --index X.gcrank --bed windows.bed [--remove-ambiguous]` that writes `bedtools nuc`-compatible columns (`pct_gc`, `num_C+num_G`, `seq_len`, ...). Also consider `--window W --step S` for genome-wide tiling. This is the main thing that would make the tool useful in pipelines. Reuse the exit-code scheme.

### P1-3. Unwrapped (single-line) FASTA costs extra memory equal to the longest line  **[verified]**
- **Where:** `src/line_reader.cpp:57`. A line that does not end inside the 64 KiB buffer is copied whole into `scratch_`.
- **Evidence:** 300 Mb genome, Release build: 481 MB peak RSS when wrapped at 60 columns, 757 MB as a single line (+276 MB, with identical output). This grows with chromosome length: a 1.5 Gb chromosome adds roughly 1.5–3 GB. Large plant genomes (*Allium*, wheat, conifers) commonly ship unwrapped.
- **Fix:** Stop materialising whole lines. Have `FastaReader` consume the buffer directly (a byte-stream state machine: at line start, `>` means header, otherwise sequence until `\n`). Only header lines need buffering. This also removes one copy for every sequence byte.

### P1-4. The parallel pipeline parallelises the cheap step; the serial producer is the bottleneck  **[verified]**
- **Evidence** (300 Mb, 2-core box, Release):

  | Config | Wall time | Peak RSS |
  |---|---|---|
  | `--threads 1` | 4.0 s | 481 MB |
  | `--threads 2` | 3.7 s | 481 MB |
  | `--threads 2 --chunk-size 10M` | 3.1 s | 138 MB |

  In-tree SHA-256 of the same file on its own: **1.4 s**, about 35% of the build, run serially in a separate pass before parsing.
- **Where the time goes:**
  1. `GcRank::build` hashes the whole file first (`gc_rank.cpp:503`), so the input is read twice. This also opens a time-of-check/time-of-use (TOCTOU) window: the file is hashed, then re-opened and parsed, so a file that changes in between gives a hash that does not match the parsed data.
  2. The producer does gzip inflate, line splitting, per-byte classification, upper-casing and `push_back` (`fasta_reader.cpp:63-87`). It then copies the whole block again into the slot (`gc_rank.cpp:154`).
  3. Workers re-classify every byte (`pack_block`). The upper-casing done in step 2 is wasted, because `base_class` is case-insensitive anyway.
- **Fixes (in order of payoff):**
  - Hash the bytes **as they are read** (wrap the `InputStream` in a hashing decorator). This is one pass and has no TOCTOU window.
  - Move validation and classification into the workers. The producer should only cut raw byte ranges at line boundaries and track contig and header events; workers validate and pack in one pass. Report the line and column of invalid bytes from the worker, which needs a per-chunk line offset.
  - `std::swap` the reader's block into the slot instead of copying it.
  - Call `gzbuffer(gz, 1<<20)` (the default is 8 KiB). Optionally support libdeflate, or parallel BGZF decompression for `.bgz` inputs, which are very common for reference genomes.
  - Lower the default `chunk_bits` from 100 M to about 8–16 M. At 100 M, a 300 Mb genome has only 3 work units and peak memory is 3.5× higher with no speed benefit.
- The README performance and memory claims should be re-measured after these changes.

### P1-5. Checksums: hand-rolled CRC-32 is ~5× slower than the zlib already linked; sidecars are never read  **[verified for speed]**
- **Evidence:** CRC-32 over 305 MB took 0.83 s in-tree versus 0.17 s with zlib's `crc32`. `crc32.h:5-7` justifies the in-tree version by "not tied to an optional dependency", but zlib is a **required** dependency (`CMakeLists.txt:27`).
- **Load cost:** `OpenOptions::validate_crc` defaults to `true` and reads every data file in full before mapping (`gc_rank.cpp:615-617`). That cancels the point of `mmap`. A ~16 Gbp genome has about 4.2 GB of index, which is roughly 12 s of CRC work on **every** load.
- **Fixes:**
  - Use zlib `crc32_z`.
  - Document the cost of `validate_crc`, or default it to `false` when `mmap = true`. Alternatively, add a cheap structural check: store `gc_total` / `valid_total` in the manifest and compare them with the last superblock rank entry (see P1-7).
  - The `*.crc32` sidecar files are written (`gc_rank.cpp:592-594`) but never read by anything. Either verify them on load as a cross-check, or drop them, since the manifest already holds the same CRCs.
- Also: `crc32_file` puts a 256 KiB `std::array` on the stack (`crc32.cpp:40`), which risks overflow on threads with small stacks (for example, macOS secondary threads get 512 KiB). Allocate it on the heap.

### P1-6. Byte-order policy is inconsistent and the cross-endian paths are untested
- `.bv` files are always little-endian, but `.rb` files are written in **host** order and tagged in the manifest (`bitvector.cpp:222`). That means two reader paths, `map_rank` and the swapped copy, and the swapped one never runs in CI.
- **Recommendation:** for format 2.0, always write `.rb` little-endian too, so readers on little-endian hosts (nearly all of them) never swap. If you keep the current format, add a big-endian CI job (for example, `s390x` via `uraimo/run-on-arch-action` or qemu-user) that loads the golden index.
- On a big-endian host, `map_bv` copies the words to the heap but still sets `mapped_ = true`, so `is_mapped()` gives the wrong answer.

### P1-7. Rank-index validation is weak when `validate_crc = false`
- `validate_rank_index` (`bitvector.cpp:340`) only checks that superblock ranks are monotonic and that basic-block entries are ≤ 4096. A corrupted `.rb` file produces silently wrong counts.
- **Fix:** store `gc_total` and `valid_total` in the manifest; check `sb.back()` against them, and check that `bb` entries are non-decreasing within each superblock. Both are O(index size / 512), which is cheap.

### P1-8. A bad `SOURCE_DATE_EPOCH` fails only after the whole build  **[verified]**
- `created_from_env()` is called at the end of `build()` (`gc_rank.cpp:529`). An invalid value wastes the entire build, which can take hours on a large genome, before failing with exit 1.
- **Fix:** validate it at the start of `build()` (and in the CLI's argument parsing).

### P1-9. Non-seekable inputs (pipes, `<(zcat …)`, `/dev/stdin`) are rejected  **[verified]**
- `fs::file_size` (`gc_rank.cpp:499`) fails with "Operation not supported" for process substitution. Once hashing is streamed (P1-4), `file_bytes` can be counted while reading, and `--fasta -` can be supported.

### P1-10. `build()` lacks the strong exception guarantee
- `build()` clears the object's state before parsing (`gc_rank.cpp:508-512`), so a failed build leaves a previously loaded index unusable. `load()` already builds into temporaries and then commits; do the same in `build()`.

### P1-11. Library writes to `stderr`
- `load_bitvector` prints a warning when it falls back from mmap (`gc_rank.cpp:476`), and verbose progress is printed directly. Libraries should not write to stderr. Add an optional log callback (`std::function<void(LogLevel, std::string_view)>`) to `BuildOptions` / `OpenOptions`, and let the CLI install the stderr sink.

---

## P2 — Maintainability, packaging, API hygiene

1. **Too much is public.** `CMakeLists.txt:91` installs every header, including `thread_pool.h`, `line_reader.h`, `input_stream.h`, `sha256.h`, `crc32.h` and `bitvector.h`. Move internals to `src/` (or `include/gcrank/detail/`, not installed). `gc_rank.h` includes `bitvector.h` → `input_stream.h`; consider a pimpl so the public header shows only `GcRank`, the options structs, `ContigInfo`, `Manifest` and `error.h`.
2. **The C++23 requirement is imposed on consumers.** `target_compile_features(gcrank PUBLIC cxx_std_23)` is stricter than needed. Once P2-1 is done, the public headers need only C++20 (`span`, `string_view`). Make it `PUBLIC cxx_std_20` and `PRIVATE cxx_std_23`. Note that the code builds with GCC 13 apart from `<print>`. Replacing `std::print` with `std::format` + `fwrite` would drop the GCC 14 requirement, which matters on HPC clusters with older toolchains.
3. **Sanitizer flags leak into the export.** `GCRANK_TSAN` adds `-fsanitize=thread` as **PUBLIC** compile and link options (`CMakeLists.txt:55-56`), so it ends up in the installed `gcrankTargets.cmake`. Make them PRIVATE, and set them on the test and CLI targets explicitly.
4. **Simplify `BuildPipeline`.** It runs N long-lived loops inside a `ThreadPool` and also keeps its own work deque, condition variable and stop flag, so it has two queues. Either submit one pool task per block, or drop `ThreadPool` and use a `std::vector<std::jthread>`. The current version has three mutexes and four condition variables for a bounded producer → workers → ordered-consumer pipeline.
5. **Inconsistent exception types:**
   - The "contig name contains a newline" check throws base `Error`, which makes the CLI exit 4 (`gc_rank.cpp:363`).
   - Querying an unbuilt index throws base `Error` (`gc_rank.cpp:635`); a `UsageError` would be more accurate.
   - `parse_contigs` runs before the CRC check in `load()` (`gc_rank.cpp:612-617`), so a corrupted `contigs.tsv` gives a confusing parse error instead of a checksum error. Check the CRC first.
6. **Defaults disagree:** `BuildOptions::threads = 4` (`gc_rank.h:33`), while the CLI uses `default_threads()` = min(8, CPUs). Use `0` to mean "auto" in the library.
7. **The reproducibility claim depends on the path string.** `manifest.source.path` stores the path exactly as given, so building the same file via a relative and an absolute path gives different bytes. (This is why the golden index has to be built from `tests/data` with a relative path.) Either store only the basename, or say this in the README's reproducibility section.
8. **JSON reader edge cases** (low risk, since the files are local):
   - Unbounded recursion on deeply nested input; add a depth limit.
   - Lone low surrogates (`\uDC00`) are accepted and emitted as invalid UTF-8.
   - `source.file_bytes` is read without a type check (`manifest.cpp:568-569`).
   - Leading zeros in numbers are accepted.
9. **CLI limits are odd:** `--threads` accepts up to 1,000,000 and `--chunk-size` up to 10¹² bases (a 1.25 TB buffer per slot). Cap them at sane values, or warn when the estimated memory exceeds available RAM.
10. **Windows paths are inconsistent:** `MappedFile` converts UTF-8 → UTF-16, but `PlainInputStream` uses `fopen` and `GzipInputStream` uses `gzopen` (use `_wfopen` / `gzopen_w`). This is low priority while Windows is officially untested.
11. **mmap hint:** call `madvise(MADV_RANDOM)` on the `.bv` mappings for random-window workloads. Document that truncating an index file while it is mapped raises `SIGBUS`. The current save path (rename) is safe; `cp` over an index is not.
12. **Nits:**
    - The `popcount64` wrapper adds nothing.
    - Use named constants for `pos >> 12` / `pos >> 9` and `static_assert` that they match `kSuperblockBits` / `kBasicBlockBits`.
    - `gc_fraction` runs `check_range` three times per call.
    - `<cstdio>` is missing from `input_stream.h`, which uses `std::FILE`.

---

## P3 — Test and CI gaps

Add tests for each P0/P1 fix, plus:

- Truncated and corrupt gzip input (P0-1); multi-member gzip; a BGZF file.
- A `ThreadPool` exception type (P0-2), if the class stays public.
- Overwrite failure that preserves the old index; concurrent saves of the same name (P0-3).
- Duplicate and empty contig names (P0-4).
- Contig-relative queries, including rejection of windows that cross contig boundaries (P1-1).
- An unwrapped FASTA with lines longer than the buffer, with a memory assertion where practical (P1-3).
- An empty FASTA with `--verbose` (P0-8).
- The documented bit layout, via an independent Python reader (P0-5).
- CI additions:
  - Clang 18 + libstdc++ (claimed but not tested).
  - macOS (claimed code path).
  - A big-endian job, or the format change in P1-6.
  - `-Werror` on the Release job.
  - Building the README example (P0-7).
  - An install + `find_package(gcrank)` consumer test, since the export is claimed but never exercised.

---

## What is already good (keep it)

- A clear on-disk contract with versioning. A golden index compared byte for byte. Reference values from an independent oracle (`data-raw/make_expected_gc.py` never calls gcrank).
- Thoughtful `bedtools nuc` parity, including the single-precision detail, with adversarial tests above 2²⁴.
- The rank structure itself (4096/512 two-level index with uint16 in-superblock ranks) is correct at every boundary I checked, compact (~4.7%), and O(1).
- A typed exception hierarchy mapped to documented CLI exit codes. Strict input validation, with helpful messages for alignment gaps.
- Sanitizer CI. Deterministic output regardless of thread count and chunk size (confirmed: single-line and wrapped inputs produce identical `.bv` files).

---

## Suggested order of work

1. P0-1, P0-2, P0-4, P0-8, each with a test. These are small and independent.
2. P0-5, P0-6, P0-7: documentation and spec corrections, plus the Python spec-checker.
3. P0-3: atomic overwrite.
4. P1-1 and P1-2: contig-relative API and `gcrank-query`. These have the most user-facing value.
5. P1-4, P1-3, P1-5: streaming hash, parser restructure, zlib CRC, smaller default chunk. Then re-benchmark and update the README numbers.
6. P2 cleanup (public header surface, CMake hygiene). Do P1-6 only if a format bump is acceptable.
7. P3 CI matrix expansion.
