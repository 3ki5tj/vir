#!/usr/bin/env python

import os, sys, glob, re, getopt
from math import *



''' extrapolate virial coefficients '''


dim = 9
nstep = -1
verbose = 0
purepy = 0
purehnc = 0
errmax = 50  # error
wreport = 0 # write a report


def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] [input.dat]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -D:     dimension
   -n:     step size of the order
   -e:     maximal number in parentheses in the error estimate
   -w:     write a report
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "D:n:M:e:wv",
        [ "dim=", "order=", "errmax=", "py", "hnc", "--report",
          "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, nstep, nmax, wreport, purepy, purehnc, verbose, errmax

  for o, a in opts:
    if o in ("-D", "--dim"):
      dim = int(a)
    elif o in ("-n", "-M", "--order"):
      nstep = int(a)
    elif o in ("-N", "--nmax"):
      nmax = int(a)
    elif o in ("-e", "--errmax"):
      errmax = int(a)
    elif o in ("-w", "--report"):
      wreport = 1
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose=",):
      verbose = int(a)
    elif o in ("--py",):
      purepy = 1
    elif o in ("--hnc",):
      purehnc = 1
    elif o in ("-h", "--help",):
      usage()

  if nstep <= 0:
    nstep = 1 if wreport else 16



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



def virextrapolate2(vir1, resol1, vir2, resol2):
  q = 1.*resol1/resol2
  dvir = (vir2 - vir1) / (q*q - 1)
  virlimit = vir1 - dvir
  return virlimit, dvir



def extrapolate(dim, order, tag = "PYc", col = 3):
  ''' estimate virial coefficients '''

  # 1. extract all virials coefficients from the sources
  template = "*Bn%sD%dn*.dat" % (tag, dim)
  fns = glob.glob(template)
  ls = []
  for fn in fns:
    m = re.search("R([0-9]*)M([0-9]*)(.*).dat", fn)
    if not m: continue
    rmax = int(m.group(1))
    npt = int(m.group(2))
    prec = m.group(3)
    try:
      arr = open(fn).readlines()[0].split()
      for x in arr:
        m = re.match(r"([0-9]+).[0-9]+", x)
        if m and int(m.group(1)) == rmax:
          rmax = float(x)
    except:
      pass
    resol = 1.*npt/rmax
    Bn = getvir(fn, order, col)
    if Bn == 0: continue
    ls += [(Bn, resol, rmax, npt, prec, fn),]

  # 2. sort the list by resolution
  ls = sorted(ls, key = lambda x: x[1])
  nls = len(ls)
  resol0, fn0, prec0, rmax0 = 0, "", "", 0
  newls = []
  for ils in range(nls):
    Bn, resol, rmax, npt, prec, fn = ls[ils]
    if fabs(resol - resol0) < 0.1:
      pdiff = precdiff(prec, prec0)
      if fn == "_" + fn0 or pdiff > 0 or (pdiff == 0 and rmax > rmax0):
        newls[-1] = (Bn, resol, rmax, npt, prec, fn)
      elif fn0 == "_" + fn or pdiff < 0 or (pdiff == 0 and rmax < rmax0):
        pass
      else:
        print "don't know how to compare %s and %s" % (fn, fn0)
        raw_input()
    else:
      newls += [(Bn, resol, rmax, npt, prec, fn),]
    resol0, prec0, fn0, rmax0 = resol, prec, fn, rmax

  ls = newls
  if len(ls) == 0:
    virlimit, err = 0, 0
  elif len(ls) == 1:
    virlimit, err = ls[0][0], 0
  elif len(ls) == 2: # two files
    vir1, resol1 = ls[-1][0], ls[-1][1]
    vir2, resol2 = ls[-2][0], ls[-2][1]
    virlimit, dvir = virextrapolate2(vir1, resol1, vir2, resol2)
    err = fabs(dvir)
  else: # three or more files
    vir1, resol1 = ls[-1][0], ls[-1][1] # the most accurate
    vir2, resol2 = ls[-2][0], ls[-2][1]
    vir3, resol3 = ls[-3][0], ls[-3][1]
    virlimit1, dvir1 = virextrapolate2(vir1, resol1, vir2, resol2)
    virlimit2, dvir2 = virextrapolate2(vir1, resol1, vir3, resol3)
    if verbose:
      print virlimit1, virlimit2, "error", fabs(virlimit1 - virlimit2), fabs(dvir1)
    err = min(fabs(virlimit1 - virlimit2), fabs(dvir1))
    virlimit = virlimit1

  if verbose:
    for x in ls:
      print x[0], x[0]-virlimit, x[1:]
  return virlimit, err, ls


def niceprint(x, err, dim = 0, n = 0, cnt = ""):
  try:
    import scifmt
    print "%4d|%4d: %30s %+20.10e %9.2e %s" % (
        dim, n, scifmt.scifmt(x, err).text(errmax = errmax),
        x, err, cnt)
  except Exception:
    print "%4d %24.14e %e %d" % (n, Bn, err, cnt)



if __name__ == "__main__":
  doargs()

  fns = ""
  tag, cols = "PYc", ((3, "self-consistent"),)
  if purepy:
    tag, cols = "PY", ((1, "compressibility"), (2, "virial"), (3, "ddP"),)
  elif purehnc:
    tag, cols = "HNC", ((1, "compressibility"), (2, "virial"),)

  nmin = nstep
  if nmin < 3: nmin = int((3 + nstep - 1)/nstep) * nstep

  datarr = []
  for col, colname in cols:
    print "column %d, %s:" % (col, colname)
    # loop over n until exhaustion
    for n in range(nmin, 10000, nstep):
      # compute the virial coefficients for dim, n
      Bn, err, ls = extrapolate(dim, n, tag, col)
      #print dim, n, Bn, err

      if len(ls) == 0: break

      datarr += [(n, col, Bn, err),]

      if not fns: # resolution + file name
        try:
          fns = "\n".join("%6s %s" % (x[1], x[-1]) for x in ls)
        except Exception:
          for x in ls: print x
          raw_input
      #print "%4d %24.14e %e %d" % (n, Bn, err, len(ls))
      niceprint(Bn, err, dim, n, len(ls))

  print fns

  # write a report
  if wreport:
    nmax = max(x[0] for x in datarr)
    src = ""
    for n in range(nmin, nmax + 1, nstep):
      svir = ""
      serr = ""
      for col, colname in cols:
        for x in datarr:
          if x[0] == n and x[1] == col:
            break
        else:
          print "n %s, col %s is missing" % (n, col)
          raw_input()
        svir += " %+22.14e" % x[2]
        serr += " %22.14e" % x[3]
      src += "%4d%s%s\n" % (n, svir, serr)
    fn = "xBn%sD%sn%s.dat" % (tag, dim, nmax)
    print "writing", fn
    open(fn, "w").write(src)
