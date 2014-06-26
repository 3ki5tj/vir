#!/bin/bash

# this function takes four parameters:
# program, dimension, order, number of evaluation points
runit () {
  # PY closure
  make $1 && ./$1 -D$2 -n$3 -M$4 --disk=2
  # HNC closure
  make $1 && ./$1 -D$2 -n$3 -M$4 --disk=2 --hnc
  # self-consistent closure
  make $1 && ./$1 -D$2 -n$3 -M$4 --disk=2 --corr
}

# compare long double and __float128
# this function takes three parameters:
# dimension, order, number of evaluation points
comp_ldbl_f128() {
  runit iegsl_l $1 $2 $3
  runit iegsl_q $1 $2 $3
}

# the following is a quick benchmark
#comp_ldbl_f128 20, 16, 512
#comp_ldbl_f128 20, 16, 768

for dim in 30 28 26 24 22 20 18 16 14 12 10; do
  comp_ldbl_f128 $dim, 128, 2048
  comp_ldbl_f128 $dim, 128, 3072
done

