#!/bin/bash

icc -DN=10 -O2 -DDG_NORING=1 -wd3180 ../../dbtool.c -o dbtool \
  && ./dbtool -i fb.bdb -j mic0/fb.bdb -k mic1/fb.bdb -o fb.bdb \
  --hash-bits=28 && \
cp fb.bdb mic0/ && \
cp fb.bdb mic1/
