#!/usr/bin/env python3
"""Compute the reference GC values for the test suite, without using gcrank.

The FASTA is scanned naively from the decompressed text, so the numbers are an
independent reference for the index queries.

Writes into tests/data/:

  expected_windows.tsv  label  start  end  gc_count  valid_count  gc_fraction  gc_fraction_remove_ambiguous
  expected_gc.tsv       scope  length  gc_count  valid_count  gc_fraction  gc_fraction_remove_ambiguous

`ERR` in a fraction column means the corresponding query must raise an error
(an empty window, or a window without A/C/G/T bases under remove_ambiguous).
Fractions are printed with repr() so an exact double comparison is possible.

Semantics implemented here (deliberately independent of gcrank): only G/C
counts as GC and only A/C/G/T counts as valid; N and the IUPAC ambiguity codes
are neither, so they stay in the default denominator but are removed by the
`remove_ambiguous` denominator. U is folded onto T (as gcrank and bedtools nuc
do), so it counts as a base but never as GC.

Usage: make_expected_gc.py [--fasta PATH] [--windows PATH] [--totals PATH]
"""
import gzip
import os
import sys

HEADER_WINDOWS = ("label", "start", "end", "gc_count", "valid_count", "gc_fraction",
                  "gc_fraction_remove_ambiguous")
HEADER_TOTALS = ("scope", "length", "gc_count", "valid_count", "gc_fraction",
                 "gc_fraction_remove_ambiguous")

# A/C/G/T, N, U and the IUPAC ambiguity codes, in either case.
ACCEPTED = set("ACGT" "NURYSWKMBDHV" "acgt" "nuryswkmbdhv")


def open_maybe_gzip(path):
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_fasta(path):
    """Returns (names, lengths, sequence) with upper-cased sequence text."""
    names = []
    lengths = []
    chunks = []
    with open_maybe_gzip(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()
                if not header:
                    raise SystemExit("empty header line in %s" % path)
                names.append(header[0])
                lengths.append(0)
                chunks.append([])
                continue
            if not chunks:
                raise SystemExit("sequence data before the first header in %s" % path)
            # Whitespace inside sequence lines is ignored by gcrank.
            sequence = "".join(line.split())
            for letter in sequence:
                if letter not in ACCEPTED:
                    raise SystemExit("invalid character %r in %s" % (letter, path))
            # Uracil is the RNA spelling of T and is folded onto it.
            sequence = sequence.upper().replace("U", "T")
            chunks[-1].append(sequence)
            lengths[-1] += len(sequence)
    return names, lengths, ["".join(chunk) for chunk in chunks]


def counts(sequence, start, end):
    """GC and valid counts; only G/C and A/C/G/T count, everything else not."""
    window = sequence[start:end]
    gc = window.count("G") + window.count("C")
    valid = window.count("A") + window.count("C") + window.count("G") + window.count("T")
    return gc, valid


def fraction(gc, denominator):
    if denominator == 0:
        return "ERR"
    return repr(gc / denominator)


def offsets(lengths):
    result = []
    running = 0
    for length in lengths:
        result.append(running)
        running += length
    return result


def build_windows(names, lengths, sequence):
    start = offsets(lengths)
    total = len(sequence)
    windows = []

    def add(label, first, last):
        gc, valid = counts(sequence, first, last)
        windows.append((label, first, last, gc, valid,
                        fraction(gc, last - first), fraction(gc, valid)))

    add("whole_genome", 0, total)
    add("first_base", 0, 1)
    add("last_base", total - 1, total)
    add("empty_window", 0, 0)
    add("superblock_0", 0, 4096)
    add("superblock_1", 4096, 8192)
    add("straddle_superblock", 4096 - 100, 4096 + 100)
    add("straddle_basic_block", 512 - 33, 512 + 33)
    add("unaligned_small", 123, 456)
    add("tail_unaligned", total - 123, total)
    add("across_first_contigs", start[1] - 50, start[2] + 37)

    for i, name in enumerate(names):
        add("contig:%s" % name, start[i], start[i] + lengths[i])

    n_index = names.index("chr_n_rich")
    add("n_rich_all_n", start[n_index], start[n_index] + 3)
    add("n_rich_mixed", start[n_index] + 3, start[n_index] + 300)
    add("n_rich_whole", start[n_index], start[n_index] + lengths[n_index])

    boundary = names.index("chr_boundary_4097")
    add("n_every_97", start[boundary], start[boundary] + 97 * 3)

    mixed = names.index("chr_mixed_case")
    add("mixed_case_half", start[mixed], start[mixed] + lengths[mixed] // 2)
    add("mixed_case_tail", start[mixed] + lengths[mixed] - 5, start[mixed] + lengths[mixed])

    iupac = names.index("chr_iupac_codes")
    plain_bases = 60
    code_bases = 120
    add("iupac_plain_part", start[iupac], start[iupac] + plain_bases)
    add("iupac_codes_part", start[iupac] + plain_bases, start[iupac] + plain_bases + code_bases)
    add("iupac_uracil_run", start[iupac] + plain_bases + code_bases,
        start[iupac] + lengths[iupac])
    add("iupac_whole", start[iupac], start[iupac] + lengths[iupac])
    add("iupac_boundary", start[iupac] + plain_bases - 10, start[iupac] + plain_bases + 10)
    return windows


def build_totals(names, lengths, sequence):
    rows = []

    def emit(scope, text):
        gc = text.count("G") + text.count("C")
        valid = text.count("A") + text.count("C") + text.count("G") + text.count("T")
        rows.append((scope, len(text), gc, valid, fraction(gc, len(text)),
                     fraction(gc, valid)))

    emit("genome", sequence)
    start = offsets(lengths)
    for i, name in enumerate(names):
        emit("contig:%s" % name, sequence[start[i]:start[i] + lengths[i]])
    return rows


def write_table(path, header, rows):
    with open(path, "w") as handle:
        handle.write("# " + "\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(field) for field in row) + "\n")
    print("wrote %s (%d rows)" % (path, len(rows)))


def main(argv):
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.normpath(os.path.join(here, os.pardir, "tests", "data"))
    fasta = os.path.join(data_dir, "test_small.fa")
    windows_path = os.path.join(data_dir, "expected_windows.tsv")
    totals_path = os.path.join(data_dir, "expected_gc.tsv")

    args = argv[1:]
    while args:
        flag = args.pop(0)
        if flag == "--fasta":
            fasta = args.pop(0)
        elif flag == "--windows":
            windows_path = args.pop(0)
        elif flag == "--totals":
            totals_path = args.pop(0)
        else:
            raise SystemExit("unknown argument: %s" % flag)

    names, lengths, chunks = parse_fasta(fasta)
    sequence = "".join(chunks)
    write_table(windows_path, HEADER_WINDOWS, build_windows(names, lengths, sequence))
    write_table(totals_path, HEADER_TOTALS, build_totals(names, lengths, sequence))
    print("scanned %d contigs, %d bases from %s" % (len(names), len(sequence), fasta))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
