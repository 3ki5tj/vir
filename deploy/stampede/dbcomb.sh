#!/bin/bash

opt=
if [ $# -ge 1 ]; then opt="-DN=$1"; fi

cc $opt -O2 ../../dbtool.c -lm && ./a.out -i fb.bdb -j mic0/fb.bdb -k mic1/fb.bdb -o fb.bdb
