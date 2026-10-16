#!/usr/bin/env python3
"""Check `pct_gc` identity between gcrank-query and `bedtools nuc`.

This is the tool test for the project's headline claim: gcrank's GC fractions
are bit-identical to `bedtools nuc`'s `pct_gc` at six printed decimals, down to
the shared single-precision arithmetic.

It builds an index for tests/data/test_small.fa, generates an exhaustive set of
BED windows (every offset for sizes 1-65, block boundaries, the ambiguity-code
and N-rich contigs, whole contigs, and random windows), and compares the two
tools row by row: `pct_gc` as printed strings, plus `num_C+num_G`,
`num_A+num_C+num_G+num_T` and `seq_len`. It also checks the adversarial
rounding cases where a double-precision quotient would print a different last
digit (and where the counts themselves exceed 2^24).

Exit codes: 0 identical, 1 mismatch or tool failure, 77 skipped (bedtools or
the CLI binaries are missing), so CTest reports a skip instead of a failure.

Usage:
  check_bedtools_parity.py [--query PATH] [--build PATH] [--bedtools PATH]
                           [--data DIR] [--tmp DIR]
"""

import argparse
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile


def run(cmd):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{proc.stderr}")
    return proc.stdout


def contigs(path):
    out = []
    name = None
    length = 0
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    out.append((name, length))
                name = line[1:].split()[0].decode()
                length = 0
            else:
                length += len(line.strip())
        if name is not None:
            out.append((name, length))
    return out


def fixture_windows(genome):
    rng = random.Random(20240607)
    windows = []
    for name, length in genome:
        for size in (1, 2, 3, 7, 63, 64, 65):
            if size <= length:
                windows.extend((name, s, s + size) for s in range(0, length - size + 1))
        for size in (100, 511, 512, 513, 4095, 4096, 4097):
            if size > length:
                continue
            starts = set(range(0, length - size + 1, 97))
            for boundary in (0, 512, 4096, 8192, 512 * 7, length - size):
                for candidate in (boundary, boundary - 1):
                    if 0 <= candidate <= length - size:
                        starts.add(candidate)
            windows.extend((name, s, s + size) for s in sorted(starts))
        windows.append((name, 0, length))
        for cut in range(1, length):
            if cut % max(1, length // 5) == 0:
                windows.append((name, 0, cut))
                windows.append((name, cut, length))
        for _ in range(3000):
            start = rng.randrange(length)
            windows.append((name, start, rng.randrange(start + 1, length + 1)))
    return windows


def parse_bedtools(text):
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        rows.append(((f[0], int(f[1]), int(f[2])), f[4],
                     int(f[6]) + int(f[7]),
                     int(f[5]) + int(f[6]) + int(f[7]) + int(f[8]), int(f[11])))
    return rows


def parse_gcrank(text):
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        rows.append(((f[0], int(f[1]), int(f[2])), f[3], int(f[4]), int(f[5]), int(f[6])))
    return rows


def compare(query, build, bedtools, fasta, index_dir, windows, label):
    bed = os.path.join(index_dir, "windows.bed")
    with open(bed, "w") as handle:
        for chrom, start, end in windows:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    theirs = parse_bedtools(run([bedtools, "nuc", "-fi", fasta, "-bed", bed]))
    ours = parse_gcrank(run([query, "--index", index_dir, "--bed", bed]))
    if len(theirs) != len(ours):
        print(f"{label}: row count differs: bedtools {len(theirs)}, gcrank {len(ours)}")
        return False
    failures = 0
    for a, b in zip(theirs, ours):
        if a[0] != b[0]:
            failures += 1
            if failures <= 5:
                print(f"{label}: interval mismatch: {a[0]} vs {b[0]}")
            continue
        for field, x, y in zip(("pct_gc", "num_gc", "num_valid", "seq_len"), a[1:], b[1:]):
            if x != y:
                failures += 1
                if failures <= 5:
                    print(f"{label}: {field} differs at {a[0]}: bedtools {x!r} vs gcrank {y!r}")
    print(f"{label}: {len(ours)} windows compared, {failures} mismatches")
    return failures == 0


def adversarial_cases(tmp):
    """Windows where float32 and double arithmetic print different decimals."""
    def f32(x):
        return struct.unpack("f", struct.pack("f", x))[0]

    def pct_float(gc, total):
        return f"{f32(f32(gc) / f32(total)):.6f}"

    def pct_double(gc, total):
        return f"{gc / total:.6f}"

    # The documented case (float prints 0.453301, double 0.453300).
    cases = [(273711, 603818)]
    # And one above 2^24, where the counts themselves round to float32; find a
    # pair whose printed values differ between the two arithmetics.
    total = (1 << 24) + 1
    for gc in range(int(total * 0.4), int(total * 0.6)):
        if pct_float(gc, total) != pct_double(gc, total):
            cases.append((gc, total))
            break
    if len(cases) == 1:
        raise RuntimeError("no adversarial case above 2^24 found")
    for i, (gc, length) in enumerate(cases):
        fasta = os.path.join(tmp, f"adv{i}.fa")
        with open(fasta, "w") as handle:
            handle.write(f">adv{i}\n")
            row = "G" * gc + "A" * (length - gc)
            for j in range(0, len(row), 60):
                handle.write(row[j:j + 60] + "\n")
        expected = pct_float(gc, length)
        assert expected != pct_double(gc, length), "case is not adversarial"
        yield f"adversarial {gc}/{length}", fasta, [(f"adv{i}", 0, length)], expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", default="build/gcrank-query")
    parser.add_argument("--build", default="build/gcrank-build")
    parser.add_argument("--bedtools", default="bedtools")
    parser.add_argument("--data", default="tests/data")
    parser.add_argument("--tmp", default=None)
    args = parser.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    query = args.query if os.path.isabs(args.query) else os.path.join(root, args.query)
    build = args.build if os.path.isabs(args.build) else os.path.join(root, args.build)
    bedtools = shutil.which(args.bedtools) or (args.bedtools if os.path.exists(args.bedtools) else None)
    if bedtools is None:
        print(f"skipped: bedtools ('{args.bedtools}') not found; install bedtools to run this check")
        return 77
    for tool in (query, build):
        if not os.path.exists(tool):
            print(f"skipped: {tool} not built; build the CLI first")
            return 77

    fasta = os.path.join(root, args.data, "test_small.fa")
    if not os.path.exists(fasta):
        print(f"skipped: fixture {fasta} not found")
        return 77

    tmp = args.tmp or tempfile.mkdtemp(prefix="gcrank-bedtools-")
    if args.tmp:
        shutil.rmtree(tmp, ignore_errors=True)  # stale indexes from a previous run
    os.makedirs(tmp, exist_ok=True)
    # Work on a copy: bedtools writes <fasta>.fai next to its input, and the
    # committed fixtures must stay untouched.
    work_fasta = os.path.join(tmp, "fixture.fa")
    shutil.copyfile(fasta, work_fasta)
    fasta = work_fasta
    ok = True
    try:
        run([build, "--fasta", fasta, "--out", tmp, "--name", "fixture"])
        index = os.path.join(tmp, "fixture.gcrank")
        ok &= compare(query, build, bedtools, fasta, index,
                      fixture_windows(contigs(fasta)), "fixture")
        for label, case_fasta, windows, expected in adversarial_cases(tmp):
            run([build, "--fasta", case_fasta, "--out", tmp,
                 "--name", os.path.basename(case_fasta)])
            case_index = os.path.join(tmp, os.path.basename(case_fasta) + ".gcrank")
            if not compare(query, build, bedtools, case_fasta, case_index, windows, label):
                ok = False
            printed = run([query, "--index", case_index, "--bed",
                           os.path.join(case_index, "windows.bed")])
            value = printed.splitlines()[1].split("\t")[3]
            if value != expected:
                print(f"{label}: gcrank printed {value}, expected {expected}")
                ok = False
    finally:
        if args.tmp is None:
            shutil.rmtree(tmp, ignore_errors=True)
    print("bedtools parity check:", "OK" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
