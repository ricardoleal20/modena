#!/bin/sh
# Clear dynamic loader/plugin/Python startup inputs before the interpreter loads.
set -eu
b07_python=$1
shift
exec /usr/bin/env -i PATH=/usr/bin:/bin:/usr/sbin:/sbin LC_ALL=C LANG=C TZ=UTC \
  RAYON_NUM_THREADS=1 RUST_TEST_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  "$b07_python" -I -S -B -X pycache_prefix=/dev/null "$@"
