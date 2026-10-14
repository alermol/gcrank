#!/bin/sh
# Regenerate the committed golden index tests/data/test_small.gcrank.
#
# The build is pinned with SOURCE_DATE_EPOCH and a relative --fasta path so
# that the test suite can rebuild the index byte-for-byte and compare it with
# the committed copy (see tests/test_roundtrip.cpp).
#
# Policy (data-raw/README.md): never regenerate the golden index just to make a
# failing test pass. Regenerate only after an intentional builder/format change.
#
# Usage: make_golden_index.sh [build-dir]   (default: <repo>/build)
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)
build_dir="${1:-$root/build}"
# Resolve a relative build directory against the caller's working directory,
# because the script changes into tests/data below.
case "$build_dir" in
  /*) ;;
  *) build_dir="$PWD/$build_dir" ;;
esac

if [ ! -x "$build_dir/gcrank-build" ]; then
  echo "error: $build_dir/gcrank-build not found; build the project first" >&2
  exit 1
fi

export SOURCE_DATE_EPOCH=1700000000
cd "$root/tests/data"
rm -rf test_small.gcrank
"$build_dir/gcrank-build" --fasta test_small.fa --out . --name test_small
echo "wrote $root/tests/data/test_small.gcrank"
