#!/usr/bin/env python



import sys, os
from math import *



def loadBn(fn, col = 1):
  ls = []
  for ln in open(fn).readlines():
    if ln.startswith("#"): continue
    a = ln.strip().split()
    ls += [(int(a[0]), float(a[col])), ]
  nmax = max(n for n, Bn in ls)
  vir = [1]*(nmax + 1)
  for n, Bn in ls: vir[n] = Bn
  return vir



def geterr(fn, fnref):
  if fnref.startswith("xBn"):
    fn, fnref = fnref, fn

  Bnref = loadBn(fnref)
  Bn = loadBn(fn) # compressibility route
  if fn.find("PYc") < 0:
    Bnv = loadBn(fn, 2) # virial route
  else:
    Bnv = None
  nmax = min(len(Bn), len(Bnref))
  s = "# %s %s\n" % (fn, fnref)
  for n in range(1, nmax):
    s1 = "%4d %24.14e %24.14e %24.14e" % (
        n, Bn[n], Bnref[n], fabs(Bn[n] - Bnref[n]) )
    if Bnv:
      s1 += " %24.14e %24.14e" % (
        Bnv[n], fabs(Bnv[n] - Bnref[n]))
    s += s1 + "\n"

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
