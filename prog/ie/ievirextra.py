#!/usr/bin/env python



''' extrapolate virial coefficients '''



import os, sys, glob, re, getopt
from math import *



dim = 9
nstep = -1
verbose = 0
errmax = 10  # error
wreport = 0 # write a report
scanall = 0 # scan all dimensions
fitorder = 3
resolmin = 200

dopy = 0
dohnc = 0
dohurst = 0
dor = 0
doir = 0
dobpgg = 0
dolamc = 0
doybg = 0
dokirk = 0

refval = 0 # reference value, for debugging


def usage():
  """ print usage and die """
  print sys.argv[0], "[Options]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -D:     dimension
   -n:     step size of the order
   -e:     maximal number in parentheses in the error estimate
   -w:     write a report
   -a:     scan all dimensions (use -aw for a report)
   -Q:     minimal resolution
   --py:   PY closure (default is the self-consistent closure)
   --hnc:  HNC closure
   --r:    the Rowlinson closure
   --ir:   the inverse Rowlinson closure
   --bpgg: the BPGG closure
   --ybg:  the YBG equation
   --kirk: the Kirkwood equation
   --lamc: density-dependent lambda
   -v:     be verbose
   -vv:    be more verbose
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "D:n:M:e:wQ:av",
        [ "dim=", "order=", "errmax=",
          "py", "hnc", "hurst", "r", "ir", "bpgg", "lamc",
          "ybg", "kirk",
          "report", "all",
          "ref=", "resmin=",
          "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, nstep, nmax, wreport, scanall, verbose, errmax
  global dopy, dohnc, dolamc, dohurst, dor, doir, dobpgg
  global doybg, dokirk
  global refval, resolmin

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
    elif o in ("-a", "--all"):
      scanall = 1
    elif o in ("-Q", "--resmin"):
      resolmin = float(a)
    elif o in ("--ref",):
      refval = float(a)
    elif o in ("--py",):
      dopy = 1
    elif o in ("--hnc",):
      dohnc = 1
    elif o in ("--lamc",):
      dolamc = 1
    elif o in ("--hurst",):
      dohurst = 1
    elif o in ("--r",):
      dor = 1
    elif o in ("--ir",):
      doir = 1
    elif o in ("--bpgg",):
      dobpgg = 1
    elif o in ("--ybg",):
      doybg = 1
    elif o in ("--kirk",):
      dokirk = 1
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  if nstep <= 0:
    nstep = 1 if wreport else 32



def getvir(fn, order, col):
  ''' the virial coefficient of `order'th order from `col'th column '''
  for s in open(fn).readlines():
    if s.strip() == "" or s[0] == '#':
      continue
    arr = s.strip().split()
    iord = int(arr[0])
    if col == -1:
      iord -= 1
    if iord == order:
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
  qq = (1.*resol1/resol2)**2
  vir12 = (vir1*qq - vir2) / (qq - 1)
  return vir12, min(fabs(vir1 - vir12), fabs(vir2 - vir12))



def linsolve(mat, b):
  ''' solve the linear equation:
            mat a = b
      by Gaussian elimination '''
  m = len(mat)
  a = b[:]
  for i in range(m):
    # choose the pivot
    im = i
    xm = fabs(mat[i][i])
    for j in range(i+1, m):
      x = fabs(mat[j][i])
      if x > xm:
        im = j
        xm = x

    # swap the row i and im
    if im != i:
      for j in range(i, m):
        mat[i][j], mat[im][j] = mat[im][j], mat[i][j]
      a[i], a[im] = a[im], a[i]
    #print mat[i], a[i]

    # normalize the ith row
    x = 1./mat[i][i]
    for j in range(i, m):
      mat[i][j] *= x
    a[i] *= x
    #print mat[i], a[i]
    #raw_input()

    # use the ith row to eliminate the kth row
    for k in range(m):
      if k == i: continue
      for j in range(i+1, m):
        mat[k][j] -= mat[k][i] * mat[i][j]
      a[k] -= mat[k][i] * a[i]
      mat[k][i] = 0;
  return a



def calcres(a, ls):
  M = len(a)
  s = res = 0
  for yx in ls:
    resol = yx[1]
    x = resol**(-2)
    y = yx[0]
    wt = resol*2
    yy = sum(x**j * a[j] for j in range(M))
    try:
      res += (yy - y)**2 * wt
    except Exception:
      print yy, y, wt
      exit()
    s += wt
  return res / s



def regression(ls, M = 3):
  ''' extrapolate the virial coefficient by regression '''
  n = len(ls)
  mat = [[0 for j in range(M)] for i in range(M)]
  b = [0 for i in range(M)]
  for a in ls:
    resol = a[1]
    x = resol**(-2) # a[1] is the resolution
    y = a[0] # the virial coefficient
    #wt = 1
    wt = resol**2
    for i in range(M):
      for j in range(M):
        mat[i][j] += x**(i+j) * wt
      b[i] += y * x**i * wt
  a = linsolve(mat, b)
  res = calcres(a, ls)
  return a[0], sqrt(res)



def searchdatalist(dim, order, tag, col):
  template = "*Bn%sD%dn*.dat" % (tag, dim)
  fns = glob.glob(template)
  ls = []
  #print fns, template
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
  return ls



def sortdatalist(ls, resolmin = 200):
  ''' sort the list by resolution '''

  # 1. sort the list by resolution, which is x[1]
  ls = sorted(ls, key = lambda x: x[1])

  # 2. remove low resolution data
  ls1 = [x for x in ls if x[1] > resolmin]
  if len(ls1) >= 4: ls = ls1

  # 3. prune data set with the same resolution
  resol0, fn0, prec0, rmax0 = 0, "", "", 0
  newls = []
  for ils in range(len(ls)):
    Bn, resol, rmax, npt, prec, fn = ls[ils]
    if fabs(resol - resol0) < resol0*0.01:
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
  nls = len(ls)
  # 4. detect bad data points that violate the monotonicity
  if nls >= 3:
    # 4.A comput the sign
    sgn = ls[-1][0] - ls[-2][0]
    if sgn > 0: sgn = 1
    else: sgn = -1

    # 4.B search the list for bad ordering
    bad = []
    for i in range(nls - 1):
      if ( (ls[i+1][0] - ls[i][0]) * sgn < 0 and
           fabs(ls[i+1][1] - ls[i][1]) > 0.01*ls[i][1] ):
        bad += [(i, i+1),]
    if bad:
      print "Warning: the list has violations: %s" % bad
      for i in range(nls): print "%4d: %s" % (i, ls[i])
      #raw_input()

  return newls



def estimatevir3(ls):
  if len(ls) == 0:
    virlimit, err = 0, 0
  elif len(ls) == 1:
    virlimit, err = ls[0][0], 0
  elif len(ls) == 2: # two files
    vir1, resol1 = ls[-1][0], ls[-1][1]
    vir2, resol2 = ls[-2][0], ls[-2][1]
    virlimit, err = virextrapolate2(vir1, resol1, vir2, resol2)
  else: # three or more files
    # 1. deterministic result from the most accurate data set
    virlimit, err0 = regression(ls[-fitorder:], fitorder)
    # 2. compute the error by fitting against all data sets
    forder = min(len(ls)-1, fitorder)
    vir, err = regression(ls, forder)
    #print forder, vir, err
    if verbose:
      print "%s %s %.10e %.10e %e %e" % (len(ls), fitorder, vir, virlimit, err, vir - virlimit)
      # print the result
      if verbose >= 2:
        # print out the results
        for x in ls:
          print x
        # print out the fitting result from different orders
        for ford in range(2, forder + 1):
          vir1, err1 = regression(ls, ford)
          print " order %4d: %.10e %.10e" % (ford, vir1, err1)

    err = max(err, fabs(virlimit - vir))
    # NOTE: here, we are using the regression result
    # but is it more reliable than the one obtained from
    # the highest-resolution data set
    virlimit = vir
  return virlimit, err



def extrapolate(dim, order, tag = "PYc", col = 3):
  ''' estimate virial coefficients '''

  # 1. extract all virials coefficients from the sources
  lsall = searchdatalist(dim, order, tag, col)

  # 2. sort the list by resolution, remove similar results
  ls0 = sortdatalist(lsall, resolmin)

  # 3. compute the virial coefficient
  virlimit0, err0 = estimatevir3(ls0)
  virlimit, err, ls, lsprec = virlimit0, err0, ls0, "all"

  # 4. classify the list by precision
  #    do the list with a particular precision
  #    x[-2] is the approximate precision
  precs = list(set(x[-2] for x in lsall))
  if len(precs) > 1:
    for prec in precs:
      if prec == "": # skip the default double precision set
        continue
      ls1 = sortdatalist([x for x in lsall if x[-2] == prec])
      if len(ls1) <= fitorder: continue
      virlimit1, err1 = estimatevir3(ls1)
      if err1 != 0 and err1 < 0.9*err:
        print "replace set [%s] by set [%s], %s(%s) -> %s(%s)" % (
            lsprec, prec, virlimit, err, virlimit1, err1)
        virlimit, err, ls, lsprec = virlimit1, err1, ls1, prec
      if verbose:
        print "list %4s %d: %24s %20s %24s %20s " % (
            prec, len(ls1), virlimit1, err1, virlimit0, err0)
  elif len(precs) == 1:
    lsprec = precs[0] # only one precision

  return virlimit, err, ls, lsprec



def niceprint(x, err, dim = 0, n = 0, cnt = "", lsprec = ""):
  try:
    import scifmt
    print "%4d|%4d: %30s %+20.10e %9.2e %s %s" % (
        dim, n, scifmt.scifmt(x, err).text(errmax = errmax),
        x, err, cnt, lsprec)
  except Exception:
    print "%4d %24.14e %e %s %s" % (n, x, err, cnt, lsprec)



def doit(dim):
  fns = ""
  tag, cols = "PYc", ((3, "self-consistent"),)
  if dopy:
    tag, cols = "PY", ((1, "compressibility"), (2, "virial"), (3, "ddP"), (-1, "cavity"))
  elif dohnc:
    tag, cols = "HNC", ((1, "compressibility"), (2, "virial"), (-1, "cavity"))
  elif dohurst:
    tag, cols = "H", ((1, "compressibility"), (2, "virial"))
  elif dor:
    tag, cols = "R", ((1, "compressibility"), (2, "virial"))
  elif doir:
    tag, cols = "IR", ((1, "compressibility"), (2, "virial"))
  elif dobpgg:
    tag, cols = "BPGG", ((1, "compressibility"), (2, "virial"))
  elif doybg:
    tag, cols = "YBG", ((1, "compressibility"), (2, "virial"), (-1, "cavity"))
  elif dokirk:
    tag, cols = "K", ((1, "compressibility"), (2, "virial"), (-1, "cavity"))
  elif dolamc:
    tag, cols = "PYl", ((3, "self-consistent"),)

  nmin = nstep
  if nmin < 3: nmin = int((3 + nstep - 1)/nstep) * nstep

  datarr = []
  for col, colname in cols:
    print "column %d, %s:" % (col, colname)
    # loop over n until exhaustion
    for n in range(nmin, 10000, nstep):
      # compute the virial coefficients for dim, n
      Bn, err, ls, lsprec = extrapolate(dim, n, tag, col)
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
      niceprint(Bn, err, dim, n, len(ls), lsprec)

  print fns

  # write a report
  if wreport and datarr:
    nmax = max(x[0] for x in datarr)
    src = ""
    for n in range(1, nmax + 1, nstep):
      svir = ""
      serr = ""
      for col, colname in cols:
        if n < nmin:
          x = [0, 0, 1.0, 0.0]
        else:
          for x in datarr:
            if x[0] == n and x[1] == col:
              break
          else:
            print "Warning: n %s, col %s is missing" % (n, col)
            x = [0, 0, 0.0, 0.0]
            #raw_input()
        svir += " %+22.14e" % x[2]
        serr += " %22.14e" % x[3]
      src += "%4d%s%s\n" % (n, svir, serr)
    fn = "xBn%sD%sn%s.dat" % (tag, dim, nmax)
    print "writing", fn
    open(fn, "w").write(src)



if __name__ == "__main__":
  doargs()
  if scanall:
    for dim in range(2, 31):
      doit(dim)
  else:
    doit(dim)


