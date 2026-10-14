# data-raw — how the test fixtures and reference values are produced

Nothing in this directory is part of the library or the CLI; it exists so that
every committed test artefact can be regenerated from scratch.

## Contents

| Script | Writes | Notes |
| --- | --- | --- |
| `make_fixtures.py` | `tests/data/test_small.fa`, `test_small_crlf.fa`, `test_small.fa.gz` | Deterministic (own LCG, fixed seed). The gzip file is written with `mtime = 0` and without a stored file name, so it is reproducible; `gzip -n -k -f` decompresses to the same bytes. |
| `make_expected_gc.py` | `tests/data/expected_windows.tsv`, `expected_gc.tsv` | Scans the FASTA naively and **never calls gcrank**, so the numbers are an independent oracle. Only `G`/`C` counts as GC and only `A`/`C`/`G`/`T` as valid; `U` is folded onto `T`. `ERR` in a fraction column means the query must raise `QueryError`. |
| `make_golden_index.sh` | `tests/data/test_small.gcrank/` | Runs the built `gcrank-build` with `SOURCE_DATE_EPOCH=1700000000` and the relative path `test_small.fa`, so `test_roundtrip.cpp` can rebuild it byte for byte. |

Typical regeneration:

```sh
python3 data-raw/make_fixtures.py
python3 data-raw/make_expected_gc.py
cmake --build build -j
sh data-raw/make_golden_index.sh build
```

## Fixture design

`test_small.fa` (11395 bases, 8 contigs) deliberately covers the boundaries of
the rank index and of the parser:

* lengths 100, 512, 513, 4096 and 4097 — basic-block and superblock edges,
* `chr_n_rich`: 1000 bases, exactly 30% `N`, including a run of `N` at offset 0
  so that an all-`N` window exists (used for the `remove_ambiguous` error case),
* `chr_boundary_4097`: an `N` every 97 bases, crossing a superblock edge,
* `chr_mixed_case`: 777 bases in mixed case (the parser upper-cases),
* `chr_iupac_codes`: 60 plain bases, then 120 bases cycling through every IUPAC
  ambiguity code (`N R Y S W K M B D H V`, half of them lower case), then 120
  uracils — so pure-ambiguous, pure-`U` and mixed windows all exist,
* header lines with descriptions and without, so the "name is the first token"
  rule is exercised,
* `test_small_crlf.fa` carries the same sequence with CRLF line endings.

`make_expected_gc.py` folds `U` onto `T` (as gcrank and `bedtools nuc` do) and
counts only `G`/`C` as GC and only `A`/`C`/`G`/`T` as valid.

## Regeneration policy for the golden index

1. **Incompatible format change** → add a new fixture (for example
   `test_small_v2.gcrank/`). Never overwrite a golden index whose layout older
   builds can no longer read, and never rewrite history for one they can.
2. **Builder bug fix** → regenerate the index (`make_golden_index.sh`), but do
   **not** regenerate the reference values: they come from an independent
   script and are the thing the builder is being checked against.
3. **Fixture content change** → regenerate the fixtures, the reference values
   and the golden index together, always in that order
   (`make_fixtures.py` → `make_expected_gc.py` → `make_golden_index.sh`). The
   reference values are still computed from the new fixture by the independent
   script, never from gcrank output.
4. **Never regenerate silently to make a failing test pass.** A byte-level
   difference between a fresh build and the golden index means either the
   builder or the golden index is wrong; decide which, on purpose, and say so
   in the commit message.

Reference values (`expected_*.tsv`) are only regenerated when the *fixtures*
change, because they are computed from the FASTA, not from gcrank.

Note for contributors: `make_golden_index.sh` pins `SOURCE_DATE_EPOCH` and uses
the relative path `test_small.fa`, because `test_roundtrip.cpp` rebuilds the
fixture with the same settings and compares the two directories byte for byte.
