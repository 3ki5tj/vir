#!/usr/bin/env python

import os, sys, glob, re, getopt
from math import *



''' extrapolate virial coefficients '''


dim = 9
nstep = 16
verbose = 0



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] [input.dat]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -D:     dimension
   -n:     step size of the order
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "D:n:M:v",
        ["dim=", "order=", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, nstep, verbose

  for o, a in opts:
    if o in ("-D", "--dim"):
      dim = int(a)
    elif o in ("-n", "-M", "--order"):
      nstep = int(a)
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose=",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()



def getvir(fn, order, col):
  for s in open(fn).readlines():
    if s.strip() == "" or s[0] == '#':
      continue
    arr = s.strip().split()
    if int(arr[0]) == order:
      return float(arr[col])
  return 0



def preccalc(prec):
  if prec == "": return 53
  if prec == "ldbl": return 64
  if prec == "f128": return 113
  if prec.startswith("p"): return int(prec[1:])
  print "unknown prec", prec
  raw_input()
  return 0



def precdiff(prec0, prec1):
  ip0 = preccalc(prec0)
  ip1 = preccalc(prec1)
  return ip0 - ip1



def extrapolate(dim, order, tag = "PYc", col = 3):
  ''' estimate virial coefficients '''
  template = "*Bn%sD%dn*.dat" % (tag, dim)
  fns = glob.glob(template)
  ls = []
  for fn in fns:
    m = re.search("R([0-9]*)M([0-9]*)(.*).dat", fn)
    if not m: continue
    rmax = int(m.group(1)) - 2
    if rmax % 2 == 1: rmax -= 1
    npt = int(m.group(2))
    prec = m.group(3)
    resol = npt/rmax
    Bn = getvir(fn, order, col)
    if Bn == 0: continue
    ls += [(Bn, resol, rmax, npt, prec, fn),]

  # sort the list by resolution
  ls = sorted(ls, key = lambda x: x[1])
  nls = len(ls)
  resol0, fn0, prec0 = 0, "", ""
  newls = []
  for ils in range(nls):
    Bn, resol, rmax, npt, prec, fn = ls[ils]
    if fabs(resol - resol0) < 0.1:
      pdiff = precdiff(prec, prec0)
      if fn == "_" + fn0 or pdiff > 0:
        newls[-1] = (Bn, resol, rmax, npt, prec, fn)
      elif fn0 == "_" + fn or pdiff < 0:
        pass
      else:
        print "don't know how to compare %s and %s" % (fn, fn0)
        raw_input()
    else:
      newls += [(Bn, resol, rmax, npt, prec, fn),]
    resol0, prec0, fn0 = resol, prec, fn

  ls = newls
  if len(ls) == 0: return 0, 0, ls
  elif len(ls) == 1: return ls[0][0], 0, ls

  vir1, resol1 = ls[-1][0], ls[-1][1]
  vir2, resol2 = ls[-2][0], ls[-2][1]
  q = 1.*resol1/resol2
  dvir = (vir2 - vir1) / (q*q - 1)
  virlimit = vir1 - dvir
  if verbose:
    for x in ls:
      print x[0], x[0]-virlimit, x[1], x[2], x[3], x[4]
  return virlimit, fabs(dvir), ls



def niceprint(x, err, dim = 0, n = 0, cnt = ""):
  try:
    import scifmt
    print "%4d|%4d: %24s %+20.10e %9.2e %s" % (
        dim, n, scifmt.scifmt(x, err).text(errmax = 10),
        x, err, cnt)
  except Exception:
    print "%4d %24.14e %e %d" % (n, Bn, err, cnt)



if __name__ == "__main__":
  doargs()
  fns = ""
  for n in range(nstep, 1000, nstep):
    Bn, err, ls = extrapolate(dim, n)
    #print dim, n, Bn, err
    if len(ls) == 0: break
    if not fns: # resolution + file name
      try:
        fns = "\n".join("%6s %s" % (x[1], x[-1]) for x in ls)
      except Exception:
        for x in ls: print x
        raw_input
    #print "%4d %24.14e %e %d" % (n, Bn, err, len(ls))
    niceprint(Bn, err, dim, n, len(ls))
  print fns
