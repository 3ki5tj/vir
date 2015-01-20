#!/bin/bash

icc -DN=14 -O2 -DDG_NORING=1 -wd3180 ../../dbtool.c -o dbtool \
  && ./dbtool -i mic0/fb.bdb -j mic1/fb.bdb -o fbmic.bdb \
  --hash-memmax=3e10 --hash-bits=24 -c0

cp fbmic.bdb mic0/fb.bdb

cp fbmic.bdb mic1/fb.bdb
