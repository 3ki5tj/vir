#!/bin/sh
export KMP_AFFINITY=compact
export OMP_NUM_THREADS=240
export MIC_ENV_PREFIX=MIC
export MIC_KMP_AFFINITY=compact
export MIC_OMP_NUM_THREADS=240

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/compilers/Intel/composer_xe_2013.5.192/compiler/lib/mic/

echo "data output to $1, random number seed $2, $LD_LIBRARY_PATH, $OMP_NUM_THREADS"

# allow 14G memory for the hash table
./a.mic -w 5 -1 1e14 -q 100000000 \
  --hash-bits=24 --hash-memmax=14e9 --hash-nocsep \
  --auto-level=9 --dbinp=fb.bdb -p3 -P $1 --rng $2

