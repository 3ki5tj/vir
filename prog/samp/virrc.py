#!/usr/bin/env python



''' compute the radius of convergence from the data collected in
    a grand canonical ensemble simulation '''



import os, sys, re, glob, getopt
import numpy as np
from math import *


dim = 0
hunt = 0
Rorder = 3
Zorder = 2
verbose = 0
nmin = -1
nmax = -1



def linsolve(a, b):
  ''' solve a x = b '''
  la = np.array(a)
  lb = np.array(b)
  sol = np.linalg.solve(la, lb)
  return list( sol )



def getrcR(x, y, ey, m = 3, name = "fb"):
  ''' fit log(y) = -c1*x - c2 - c3/x - ...'''
  n = len(x)
  lny = [ log(abs(yi)) for yi in y ]
  elny = [ ey[i]/(abs(y[i])) for i in range(n) ]
  mat = [None] * m
  for i in range(m): mat[i] = [0] * m
  vec = [0] * m
  for i in range(n):
    xi = x[i]
    w = 1./(elny[i] * elny[i])
    for u in range(m):
      for v in range(m):
        mat[u][v] += pow(xi, 2. - u - v) * w
      vec[u] += lny[i] * pow(xi, 1 - u) * w
  c = linsolve(mat, vec)
  s = "%14.7e*n" % c[0]
  if m >= 2: s += " %+14.7e" % c[1]
  if m >= 3: s += " %+14.7e/n" % c[2]
  for i in range(3, m): s += " %+14.7e/n**%d" % (c[i], i - 1)
  print "log|%s| = %s" % (name, s)
  return exp(-c[0])



def getZa(x, y, ey, m = 2):
  ''' y = (a1 * x + a2 + ...) / (1 - b1 / x - ...)
      y = a1 * x + a2 + ... + b1 y / x + ... '''
  mm = 2*m - 1 # total number of parameters
  mat = [None] * mm
  for i in range(mm): mat[i] = [0] * mm
  vec = [0] * mm
  n = len(x)
  for i in range(n):
    xi = x[i]
    yi = y[i]
    if ey[i] <= 0:
      print "i %d n %d, ey is zero y %g, ey %g" % (i, xi, yi, ey[i])
      continue
    w = 1/(ey[i]*ey[i])
    for u in range(m):
      for v in range(m):
        mat[u][v] += pow(xi, 2. - u - v) * w
      vec[u] += yi * pow(xi, 1 - u) * w
      for v in range(m-1):
        mat[u][v + m] += yi * pow(xi, - u - v) * w
        mat[v + m][u] += yi * pow(xi, - u - v) * w
    for u in range(m - 1):
      for v in range(m - 1):
        mat[u + m][v + m] += yi * yi * pow(xi, - 2 - u - v) * w
      vec[u + m] += yi * yi * pow(xi, -1 - u) * w
  a = linsolve(mat, vec)
  s = "%14.7e*n" % a[0]
  if m >= 2: s += " %+14.7e" % a[1]
  if m >= 3: s += " %+14.7e/n" % a[2]
  for i in range(3, m): s += " %+14.7e/n**%d" % (a[i], i - 1)
  t = "1"
  for i in range(m - 1): t += " %+14.7e/n**%d" % (-a[i + m], i + 1)
  print "Zr(n) = (%s)/(%s)" % (s, t)
  return a[0]



def getRc(narr, fbarr, fberr, fbnrarr, fbnrerr, Zrarr, eZrarr,
          Rorder = 3, Zorder = 2, hist = 0):
  ''' estimate the radius of convergence from the ring route
      and from the partition function (activity) route '''
  if fbarr != None and fbnrarr != None:
    rcR = getrcR(narr, fbnrarr, fbnrerr, Rorder, "fb/nr")
    rcZ = getrcR(narr, fbarr, fberr, Rorder, "fb")
  else:
    rcR = rcZ = 1
  b = log(rcZ)
  if Zrarr != None:
    a = getZa(narr, Zrarr, eZrarr, Zorder)
  else:
    a = 1
  rcZ /= a
  print "%s lines, %9.2e per vircoef. rcR %s(%s), rcZ %s(%s, %s)" % (
      len(narr), hist, rcR, log(rcR), rcZ, a, b)
  return rcR, rcZ



def getcol(lines, c):
  ''' get a column '''
  return [d[c] for d in lines]



def getrcZr(fn):
  ''' estimate the radius of convergence from file fn '''
  if fn.endswith(".data"):
    # try to call sum.py
    try:
      import virsum
      virsum.aggregate(None)
    except ImportError:
      pass

  isZrh = fn.startswith("Zrh")
  width = 19
  if isZrh: width = 23

  lines = open(fn).readlines()[1:]
  if nmin >= len(lines) - 3:
    print "only %d lines, bad nmin %d, nmax %d" % (
        len(lines), nmin, nmax)
    raise Exception
  lines = lines[nmin - 1: nmax]
  n = len(lines)
  # parse each line
  for i in range(n):
    lines[i] = [float(x) for x in lines[i].split()]
  narr = getcol(lines, 0)
  Zrarr = getcol(lines, 1)
  if isZrh:
    hist = getcol(lines, 5)
    nrarr = getcol(lines, 13)
    fbarr = getcol(lines, 20)
    BRarr = getcol(lines, 22)
  else:
    hist = getcol(lines, 3)
    nrarr = getcol(lines, 9)
    fbarr = getcol(lines, 16)
    BRarr = getcol(lines, 18)
  avhist = sum(hist)/n
  if len(lines[1]) > width:
    BRerr = getcol(lines, 19)
    Zrerr = getcol(lines, 21)
  else:
    BRerr = [1] * n
    Zrerr = [1] * n
  fberr = [0]*n
  fbnrarr = [0]*n
  fbnrerr = [0]*n  # fb/nr
  nrhaszero = 0
  for i in range(n):
    if nrarr[i] != 0:
      fbnrarr[i] = fbarr[i]/nrarr[i]
      fberr[i] = BRerr[i]/BRarr[i]*fbarr[i]
      fbnrerr[i] = fberr[i]/nrarr[i]
    else:
      nrhaszero = 1
      fbnrarr[i] = 0
      fberr[i] = 1
      fbnrerr[i] = 1

  if nrhaszero:
    fbarr = fbnrarr = None

  if verbose:
    for i in range(n):
      print "%3d: vir %+20.10e (%9.2e)" % (
          narr[i], BRarr[i], BRerr[i])

  return getRc(narr, fbarr, fberr, fbnrarr, fbnrerr, Zrarr, Zrerr,
        Rorder, Zorder, hist = avhist)



def huntmr(d, n, root = "."):
  """ hunt for mr files for order n """

  try:
    # try to use virsum because it produces error estimates
    import virsum
    virsum.verbose = verbose
    (x, err, tot, strtot) = virsum.dodirs(None, n, sum3 = 1)
    if x == None: return None, None, None, None, None
    else: return x[1], err[1], x[2], err[2], tot
  except ImportError:
    print "cannot import virsum.py, use the default method"
    pat = re.compile(r"mrD%dn%d.dat[0-9]+" % (d, n))
    ls = []
    for r, ds, fs in os.walk(root):
      ls += [os.path.join(r, f) for f in fs if pat.match(f)]
    fbcnt = fbsum = 0.0
    nrcnt = nrsum = 1e-30
    if len(ls) == 0:
      return None, None, None, None
    for f in ls:
      arr = open(f).readlines()[1].split()
      fbcnt += float(arr[0])
      fbsum += float(arr[1]) * float(arr[0])
      nrcnt += float(arr[2])
      nrsum += float(arr[3]) * float(arr[2])
    fb = fbsum/fbcnt
    nr = nrsum/nrcnt
    return fb, 1, nr, 1, fbcnt



def huntrc(d, nmin, nmax):
  n = max(nmin, 4)
  narr = []
  fbarr = []
  fberr = []
  fbnrarr = []
  fbnrerr = []
  totarr = []
  # information line of the output file
  src = "#  n         <fb>              e(fb)           <nr>             e(nr)             tot\n"
  while 1:
    fb, errfb, nr, errnr, tot = huntmr(d, n)
    if fb == None: break
    narr += [n,]
    fbarr += [fb,]
    fbnrarr += [fb/nr,]
    if errfb == 1 and errnr == 1:
      fberr += [1,]
      fbnrerr += [1,]
    else:
      fberr += [errfb,]
      fbnrerr += [errfb/nr,]
    totarr += [tot,]
    line = "%4d %+22.14e %9.2e %22.14e %9.2e %20.0f\n" % (
        n, fb, errfb, nr, errnr, tot)
    src += line
    print line,
    if verbose: raw_input("press Enter to continue")
    n += 1
    if n > nmax: break
  fnout = "fbnr.dat"
  open(fnout, "w").write(src)
  print "writing %s" % fnout
  return getRc(narr, fbarr, fberr, fbnrarr, fbnrerr, None, None,
        Rorder, Zorder, hist = 0)



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] [input.dat]"
  print """
  Compute the radius of convergence from a grand canonical ensemble simulation
    mcgc2.c, mcgcr2.c
  The input are something like
    ZrhD8n64.dat, ZrrD100r4n64.dat
  Without the input, the program will try search proper input dat files

  OPTIONS:
   -m: search for mrDXnX.datnnn for fb and nr
   -d: dimension
   -n: nmin
   -N: nmax
   -o: order (default is 2)
   -v: verbose
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "vmd:n:N:o:O:",
         ["hunt", "dim=", "nmin=", "nmax=",
           "order=", "Order=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, nmin, nmax, verbose, Rorder, Zorder, hunt

  fn = None

  for o, a in opts:
    if o in ("-m", "--hunt"):
      hunt = 1
    elif o in ("-d", "--dim"):
      dim = int(a)
    elif o in ("-n", "--nmin"):
      nmin = int(a)
    elif o in ("-N", "--nmax"):
      nmax = int(a)
    elif o in ("-O", "--Order"):
      Rorder = int(a)
    elif o in ("-o", "--order"):
      Zorder = int(a)
    elif o in ("-v",):
      verbose += 1
    elif o in ("-h", "--help",):
      usage()

  # check if to turn on the hunting mode
  if hunt == 0:
    curdir = os.getcwd()
    if re.search(r"D[0-9]+$", curdir):
      hunt = 1

  if hunt: # hunting mode
    if dim == 0: # guess the dimension
      m = re.match(".*D([0-9]+)", os.getcwd())
      if not m:
        print "cannot determine the dimension"
        usage()
      dim = int(m.group(1))
    print "hunting mode for D = %s" % dim
    if nmin <= 0: nmin = 4
    if nmax <= 0: nmax = 64
  else: # normal mode
    if len(args):
      fn = args[0]
    else:
      fn = glob.glob("ZrD*.data")
      if len(fn) == 0:
        fn = glob.glob("ZrhD*.data")
        if len(fn) == 0:
          print "cannot find data file"
          usage()
      fn = fn[0]
    if nmin <= 0: nmin = 20
    if nmax <= 0: nmax = 64

  return fn



if __name__ == "__main__":
  fn = doargs()
  if hunt:
    huntrc(dim, nmin, nmax)
  else:
    getrcZr(fn)
