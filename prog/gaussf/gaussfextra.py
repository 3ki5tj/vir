#!/usr/bin/env python


''' extrapolate over the number of edges to get the best estimated virial coefficient '''



import glob, os, sys
from math import *
import scifmt



def getvir(fn, dim):
  for ln in open(fn).readlines():
    s = ln.strip()
    if s.startswith("#"): continue
    d, vir = s.split()
    if int(d) == dim:
      return float(vir)
  print "cannot find D %d dimensional value from %s" % (dim, fn)
  return 0



def extrapolate(dim, n):
  # 1. find files
  fnls0 = glob.glob("gaussfn%de*E*ldblmpf.dat" % (n))
  fnls = []
  for ed in range(13, 100):
    fn = "gaussfn%de%dE%dldblmpf.dat" % (n, ed, ed)
    if os.path.exists(fn):
      try:
        fnls0.remove(fn)
      except Exception:
        print fn, fnls0
      fnls += [(ed, fn),]
  fnls = [(0, fn) for fn in fnls0] + fnls
  #print fnls

  # 2. collect virial coefficients
  ls = []
  x = 0
  for ed, fn in fnls:
    vir = getvir(fn, dim)
    x += vir
    ls += [(ed, vir, x),]
  #print ls

  # 3. extrapolate using a geometric series
  Bnprev = prev = 0
  for ed, vir, svir in ls:
    if ed == 0: continue
    rat = prev / vir
    dx = vir/(rat - 1)
    if rat != 0:
      Bn = svir + dx
    else:
      Bn = svir
    if Bnprev > 0:
      err = fabs(Bn - Bnprev)
    else:
      err = fabs(dx)
    if prev != 0:
      print "D %2d, n %2d, e %2d %-28s %.20f %s | %s rat %s" % (
          dim, n, ed, scifmt.scifmt(Bn, err).text(errmax = 100), Bn, err, svir, rat)
    prev = vir
    Bnprev = Bn

  try:
    print "Final\nD %2d, n %d, %-28s %.20f %s" % (
        dim, n, scifmt.scifmt(Bn, err).text(errmax = 100), Bn, err)
  except Exception:
    print "D %d, n %d, %.20f %s" % (dim, n, Bn, err)



if __name__ == "__main__":
  dim = 8
  n = 14
  if len(sys.argv) >= 3:
    dim = int(sys.argv[1])
    n = int(sys.argv[2])
  extrapolate(dim, n)
