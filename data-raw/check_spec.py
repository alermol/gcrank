#!/usr/bin/env python3
"""Decode a .gcrank index using only the documented bit layout.

This is a doc-test for the on-disk format. It reads genome.bv, valid.bv and
contigs.tsv with a from-scratch reader that implements exactly the rule written
in include/gcrank/format.h (and README.md):

    the .bv file is a sequence of 64-bit little-endian words; bit `i` of a
    bitvector is bit `63 - i % 64` of word `i / 64`

and checks the decoded counts against tests/data/expected_gc.tsv, whose values
are produced by make_expected_gc.py from the FASTA and never by gcrank. If the
documentation and the bytes ever disagree again, this script fails.

Usage: check_spec.py [index_dir] [expected_gc.tsv]
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INDEX = ROOT / "tests" / "data" / "test_small.gcrank"
DEFAULT_EXPECTED = ROOT / "tests" / "data" / "expected_gc.tsv"


def bitvector_words(path):
    """The .bv bytes as the documented little-endian 64-bit words."""
    data = path.read_bytes()
    if len(data) % 8 != 0:
        raise SystemExit(f"{path}: size {len(data)} is not a whole number of words")
    return [int.from_bytes(data[i:i + 8], "little") for i in range(0, len(data), 8)]


def bit(words, i):
    """Bit `i` of the vector, per the documented rule (0 = first bit)."""
    return (words[i // 64] >> (63 - (i % 64))) & 1


def bits_set(words, start, end):
    return sum(bit(words, i) for i in range(start, end))


def load_contigs(path):
    contigs = []
    for line in path.read_text().splitlines():
        if not line:
            continue
        name, length, offset = line.rsplit("\t", 2)
        contigs.append((name, int(length), int(offset)))
    return contigs


def load_expected(path):
    rows = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        rows[fields[0]] = (int(fields[1]), int(fields[2]), int(fields[3]))
    return rows


def main():
    index = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INDEX
    expected = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_EXPECTED

    genome = bitvector_words(index / "genome.bv")
    valid = bitvector_words(index / "valid.bv")
    contigs = load_contigs(index / "contigs.tsv")
    rows = load_expected(expected)

    failures = []

    def check(scope, length, gc, valid_count):
        if scope not in rows:
            failures.append(f"{scope}: missing from {expected.name}")
            return
        want = rows[scope]
        got = (length, gc, valid_count)
        if got != want:
            failures.append(f"{scope}: decoded {got}, expected {want}")

    running_offset = 0
    for name, length, offset in contigs:
        if offset != running_offset:
            failures.append(f"{name}: offset {offset}, expected {running_offset}")
        running_offset += length
    total_length = running_offset

    check("genome", total_length, bits_set(genome, 0, total_length),
          bits_set(valid, 0, total_length))
    for name, length, offset in contigs:
        check(f"contig:{name}", length, bits_set(genome, offset, offset + length),
              bits_set(valid, offset, offset + length))

    if failures:
        for problem in failures:
            print(f"spec check FAILED: {problem}", file=sys.stderr)
        raise SystemExit(1)
    print(f"spec check OK: {total_length} bits decoded per the documented rule "
          f"({len(contigs)} contigs)")


if __name__ == "__main__":
    main()
