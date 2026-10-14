#!/usr/bin/env python3
"""Generate the deterministic FASTA fixtures used by the gcrank test suite.

Writes into tests/data/:

  test_small.fa       LF line endings, ~11 kB
  test_small_crlf.fa  byte-identical sequence with CRLF line endings
  test_small.fa.gz    gzip of test_small.fa with no timestamp and no stored
                      file name (reproducible: `gzip -n -k -f` gives the same
                      decompressed bytes, though not necessarily the same
                      compressed bytes as Python's zlib)

The sequence comes from a small deterministic LCG, so regenerating the
fixtures on any machine produces byte-identical files.

Contigs cover the rank-index boundaries (100, 512, 513, 4096, 4097 bases), one
contig holds exactly 30% N, one mixes upper and lower case, one carries an N
every 97 bases so that N handling can be checked across block boundaries, and
one holds every IUPAC ambiguity code followed by a run of uracils.
"""
import gzip
import os
import sys

MASK64 = (1 << 64) - 1
SEED = 20231114

BASES = "ACGT"
# The IUPAC ambiguity codes (N, U and the ambiguity codes are handled
# separately: U is folded onto T, N is in AMBIGUOUS).
CODES = "NRYSWKMBDHV"


class Lcg:
    """Deterministic pseudo random generator (Knuth LCG)."""

    def __init__(self, seed):
        self.state = seed & MASK64

    def next(self):
        self.state = (6364136223846793005 * self.state + 1442695040888963407) & MASK64
        return self.state

    def base(self):
        return BASES[(self.next() >> 33) % 4]


def make_contigs():
    """Returns a list of (name, sequence) with deterministic content."""
    rng = Lcg(SEED)
    contigs = []

    def plain(name, length, n_every=None):
        out = []
        for i in range(length):
            if n_every is not None and i % n_every == 0:
                out.append("N")
            else:
                out.append(rng.base())
        return name, "".join(out)

    def n_rich(name, length):
        # Exactly 30% N: positions i with i % 10 < 3.
        out = []
        for i in range(length):
            out.append("N" if (i % 10) < 3 else rng.base())
        return name, "".join(out)

    def mixed_case(name, length):
        out = []
        for i in range(length):
            base = rng.base()
            out.append(base.lower() if (i // 100) % 2 == 0 else base)
        return name, "".join(out)

    def iupac(name, plain_bases, code_bases, uracil_bases):
        # Plain DNA first (so windows without any code exist), then every
        # ambiguity code (half of it lower case), then a run of uracils.
        out = [rng.base() for _ in range(plain_bases)]
        for i in range(code_bases):
            code = CODES[i % len(CODES)]
            out.append(code.lower() if (i // 55) % 2 == 0 else code)
        for i in range(uracil_bases):
            out.append("u" if i % 2 == 0 else "U")
        return name, "".join(out)

    contigs.append(plain("chr_boundary_100", 100))
    contigs.append(plain("chr_boundary_512", 512))
    contigs.append(plain("chr_boundary_513", 513))
    contigs.append(plain("chr_boundary_4096", 4096))
    contigs.append(plain("chr_boundary_4097", 4097, n_every=97))
    contigs.append(n_rich("chr_n_rich", 1000))
    contigs.append(mixed_case("chr_mixed_case", 777))
    contigs.append(iupac("chr_iupac_codes", 60, 120, 120))
    return contigs


def render(contigs, newline, wrap):
    lines = []
    for i, (name, sequence) in enumerate(contigs):
        lines.append(">%s description of %s" % (name, name) if i % 2 == 0 else ">%s" % name)
        for start in range(0, len(sequence), wrap):
            lines.append(sequence[start:start + wrap])
    return newline.join(lines) + newline


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(here, os.pardir, "tests", "data")
    data_dir = os.path.normpath(data_dir)

    contigs = make_contigs()
    lf = render(contigs, "\n", 60)
    crlf = render(contigs, "\r\n", 60)

    targets = {
        "test_small.fa": lf.encode("ascii"),
        "test_small_crlf.fa": crlf.encode("ascii"),
    }
    for name, payload in sorted(targets.items()):
        path = os.path.join(data_dir, name)
        with open(path, "wb") as handle:
            handle.write(payload)
        print("wrote %s (%d bytes)" % (path, len(payload)))

    # mtime=0 strips the timestamp; Python's gzip stores no file name.
    gz_path = os.path.join(data_dir, "test_small.fa.gz")
    with open(gz_path, "wb") as handle:
        handle.write(gzip.compress(lf.encode("ascii"), compresslevel=6, mtime=0))
    print("wrote %s (%d bytes)" % (gz_path, os.path.getsize(gz_path)))

    total = sum(len(sequence) for _, sequence in contigs)
    print("contigs: %d, bases: %d" % (len(contigs), total))
    return 0


if __name__ == "__main__":
    sys.exit(main())
