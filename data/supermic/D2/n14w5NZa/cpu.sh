#!/bin/sh

export OMP_NUM_THREADS=20
export KMP_AFFINITY=verbose,scatter

echo "random number seed $1"

export dir1=`pwd`
echo "cpu.sh: $dir1, seed $1" >> cur.dir

# allow 60 gigabytes for the hash table
./a.out -w 5 -1 1e14 -q 1000000000 --rng $1 \
  --hash-mode=1 --hash-bits=28 --hash-memmax=6e10 \
  --hash-nocsep --auto-level=9 --dbez -p3 &

