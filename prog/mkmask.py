#!/usr/bin/env python

''' generate bits mask '''

def bitsmask(n):
  return (1 << n) - 1


def printbitsmask(nmax, inv = 0):
  n = 1
  mask = bitsmask(nmax)
  for i in range(nmax / 4):
    ii = i * 4
    line = ""
    for j in range(4):
      ij = ii + j + 1
      x = bitsmask(ij)
      if inv: x = x ^ mask
      if nmax == 32:
        sval = "0x%08x" % x
      else:
        sval = "CU64(" + ("0x%016x" % x) + ")"
      line += " " + sval + ","
    print line
  print ""

printbitsmask(32, inv = 0)
printbitsmask(32, inv = 1)
printbitsmask(64, inv = 0)
printbitsmask(64, inv = 1)

