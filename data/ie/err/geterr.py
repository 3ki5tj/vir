#!/usr/bin/env python



import sys, os
from math import *



def loadBn(fn):
  ls = []
  for ln in open(fn).readlines():
    if ln.startswith("#"): continue
    a = ln.strip().split()
    ls += [(int(a[0]), float(a[1])), ]
  nmax = max(n for n, Bn in ls)
  vir = [1]*(nmax + 1)
  for n, Bn in ls: vir[n] = Bn
  return vir



def geterr(fn, fnref):
  if fnref.startswith("xBn"):
    fn, fnref = fnref, fn

  Bnref = loadBn(fnref)
  Bn = loadBn(fn)
  nmax = min(len(Bn), len(Bnref))
  s = "# %s %s\n" % (fn, fnref)
  for n in range(1, nmax):
    s += "%4d %24.14e %24.14e %24.14e\n" % (
        n, Bn[n], Bnref[n], fabs(Bn[n] - Bnref[n]))

  # output name
  fnout = fn
  if fn.endswith(".dat"): fnout = fn[:-4] + ".err"
  else: fnout = fn + ".err"
  print "saving data to %s" % fnout
  open(fnout, "w").write(s)



if __name__ == "__main__":
  fn = "xBnPYcD15n128.dat"
  fnref = "BnD15n64.dat"
  if len(sys.argv) >= 3:
    fn, fnref = sys.argv[1], sys.argv[2]
  geterr(fn, fnref)
