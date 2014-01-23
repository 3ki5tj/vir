#!/usr/bin/env python


'''
    summarize multiple thread/core data generated by
    hsrh.c, mcrat0.c, mcgc2.c, and mcgcr2.c
    X.dat, X.dat1, X.dat2, ...  --> X.data

    program     prefix    input(e.g.)       python class
    ---------------------------------------------------
    hsrh.c        int     intD3n6.dat         INT
    mcrat0.c      mr      mrD3n10.dat         MR
    mcgc2.c       Zrh     ZrhD10n64.dat       GC
    mcgcr2.c      Zrr     ZrrD100r4n64.dat    GC
'''


import os, sys, re, glob, getopt
from math import *



verbose = 1 # verbose level
rmbad = 0 # detect and remove bad data that significantly deviate from the average
mulsig = 3 # multiple of standard deviation for the above exclusion
dirnamecheck = True



def errarr(arr):
  ''' calculate the error of the array `arr' '''
  m = len(arr) # number of copies
  if m <= 1:
    print "only %d copy: %s" % (m, arr)
    raise Exception
  n = len(arr[0])
  err = [0] * n
  for i in range(n): # loop over arrays
    smx = smx2 = 0
    xb = arr[0][i]
    for k in range(m): # loop over copies
      x = arr[k][i] - xb
      smx += x
      smx2 += x * x
    avx = smx / m
    smx2 = max(smx2 / m - avx * avx, 0)
    err[i] = sqrt(smx2 / (m - 1))
  return err



def varanalysis(arr, fns = None, name = "", nsig = 3, nave = 0, verbose = True):
  ''' variance analysis: list large deviation data directories '''
  n = len(arr)
  ave = 1. * sum(arr) / n
  var = sum([x*x for x in arr]) / n - ave * ave
  if var < 0: var = 0
  sig = sqrt(var)
  badls = []
  threshold = max(fabs(ave) * nave, nsig * sig)
  for i in range(n):
    if fabs(arr[i] - ave) > threshold:
      badls += [i,]
  if len(badls) and verbose:
    if fns:
      w = max([len(fn) for fn in fns]) + 4
    else:
      w = 10
      fns = [ "data%d" % i for i in range(n) ]
    for i in badls:
      print "%simbalance: %-*s %+10.3e, %+10.3e/%7.1e (ave/sig)" % (
          name + " ", w, fns[i], arr[i], ave, sig)
  return ave, sig, badls





class INT:
  ''' a class for combining data from hsrh.c '''
  once = 0

  def __init__(me, fn):
    s = open(fn).readlines()
    me.initdata(s[0])
    me.loaddata(s[1:])



  def initdata(me, s):
    arr = s.split()
    me.tag = arr[0][1:].strip()
    me.D = int(arr[1])
    me.n = int(arr[2])
    me.m = int(arr[3])
    me.acc = [0] * me.m
    me.fac = [0] * me.m
    me.names  = [""] * me.m
    me.vir = 0
    me.evir = 0



  def loaddata(me, s):
    me.tot = float(s[0])
    for i in range(me.m):
      x = s[i + 1].split()
      me.acc[i] = float(x[0])
      me.fac[i] = float(x[1])
      me.names[i] = x[2]



  def absorb(a, b):
    ''' a += b '''
    a.tot += b.tot
    for i in range(a.m):
      a.acc[i] += b.acc[i]



  def recompute(me):
    ''' recompute vir '''
    me.vir = 0
    for i in range(me.m):
      me.vir += me.acc[i] * me.fac[i]
    me.vir /= me.tot
    for i in range(1, me.n + 1):
      me.vir /= i
    me.vir *= -(me.n - 1) * pow(2, me.n - 1);



  def save(me, fn):
    s = "#0 %d %d %d V0\n" % (me.D, me.n, me.m)
    s += "%s\n" % me.tot
    for i in range(me.m):
      s += "%16.0f %+.14e %-6s\n" % (
          me.acc[i], me.fac[i], me.names[i])
    s += "%.14e %.1e\n" % (me.vir, me.evir)
    open(fn, "w").write(s)



  @staticmethod
  def sumdat(fnbas, fnls):
    ''' sum over data in fnbas + fnls for int*.dat '''
    intg = INT(fnbas)
    nls = len(fnls) + 1
    intgls = [None] * nls
    fnls = [fnbas,] + fnls
    # sum over all file lists
    for i in range(nls):
      fn = fnls[i]
      intgls[i] = INT(fn)
      intgls[i].recompute() # compute individual vir
      if i > 0: intg.absorb(intgls[i])

    # compute the global data
    intg.recompute()

    # compute errors
    virls = [None] * nls
    for i in range(nls): # loop over copies
      virls[i] = [intgls[i].vir,]
    intg.evir = errarr( virls )[0]

    fnout = fnbas + "a"
    intg.save(fnout)
    print "saving %s, vir %.8e, %+9.2e, tot %g" % (
        fnout, intg.vir, intg.evir, intg.tot)
    return ([intg.vir,], [intg.evir,], intg.tot)





def toarr(b):
  a = [float(x) for x in b]
  if a[0] <= 0: a[0] = 1e-20
  for i in range(1, len(a)):
    a[i] *= a[0]
  return a



def xinc(a, b):
  if len(a) != len(b):
    print "arrays mismatch\n%s\n%s" % (a, b)
    raise Exception
  for i in range(len(b)):
    a[i] += b[i]
  return a



def getver(sver):
  if sver.startswith("V"):
    return int(sver[1:])
  else:
    return int(sver)



def findfile(fn):
  ''' find fn '''
  # current directory
  fn0 = fn
  if os.path.exists(fn): return fn
  # parent directory
  fn = os.path.join(os.pardir, fn0)
  if os.path.exists(fn): return fn
  # directory of the script
  d = os.path.split(os.path.realpath(__file__))[0]
  fn = os.path.join(d, fn0)
  if os.path.exists(fn): return fn
  print "Error: cannot find %s, src dir %s, %s" % (fn, d, __file__)
  return None



def findBring():
  ''' find Bring.dat '''
  return findfile("Bring.dat")



def findZn():
  ''' find Z.dat '''
  return findfile("Z.dat")



class MR:
  ''' a class for combining data from mcrat.c '''
  once = 0

  def __init__(me, fn, fnBring = None, use3 = 0):
    me.fn = fn
    s = open(fn).readlines()
    if len(s) <= 1:
      print "file %s is corrupted, delete it!\n" % os.path.realpath(fn)
      raise Exception
    me.initdata(s[0])
    me.loaddata(s[1])
    me.loadBring(fnBring)
    me.use3 = use3



  def loadBring(me, fn):
    ''' load virial coefficients from Bring.dat '''
    if not fn: fn = findBring()
    if not fn:
      print "cannot find Bring.dat"
      return

    ls = open(fn).readlines()
    for s in ls:
      bs = s.split()
      if int(bs[0]) == me.D:
        break
    else:
      print "%s: no data for D %d" % (fn, me.D)
      return
    me.Bring = fabs( float(bs[ me.n ]) )
    if not GC.once:
      if verbose:
        print "loaded n %d, D %d from %s, %s" % (me.n, me.D, fn, me.Bring)
      GC.once = 1


  Znsrc = None
  def loadZn(me, fn = None):
    ''' override the partition function if possible '''
    # the precision of the data below is about 1%
    if MR.Znsrc == None:
      if not fn: fn = findZn()
      if not fn:
        print "cannot find Zn.dat"
        return
      MR.Znsrc = open(fn).readlines()
    for s in MR.Znsrc:
      bs = s.split()
      if int(bs[0]) == me.D and len(bs) > me.n:
        me.Zn = float(bs[ me.n ])
        break
    else:
      pass
      #print "%s: no data for D %d, n %d" % (fn, me.D, me.n)



  def initdata(me, s):
    arr = s.split()
    me.tag = arr[0][1:].strip()
    me.D = int(arr[1])
    me.n = int(arr[2])
    me.version = int(arr[3][1:])
    if me.tag == 'M':
      me.Z = [1, float(arr[4]), float(arr[5])]
    elif me.version >= 1:
      me.Zn = float(arr[5])
      me.loadZn() # override Zn if possible
    me.fbsm = [[0,0], [0,0], [0,0]]
    me.nzsm = [0,0]  # only for me.tag == 'M'
    me.nrsm = [0,0]
    me.tacc = [[0,0], [0,0]]
    me.evir = None
    me.tot = 0



  def loaddata(me, s):
    x = s.strip().split()
    me.fbsm[0] = toarr((x[0], x[1]))
    me.tot = me.fbsm[0][0]
    if me.tag == 'M':
      me.fbsm[1] = toarr((x[2], x[3]))
      me.fbsm[2] = toarr((x[4], x[5]))
      me.nzsm = toarr((x[6], x[7]))
      me.nrsm = toarr((x[8], x[9]))
      me.tacc[0] = toarr((x[10], x[11]))
      me.tacc[1] = toarr((x[12], x[13]))
    else:
      me.nrsm = toarr((x[2], x[3]))



  def recompute(me):
    ''' recompute virial '''
    if me.version >= 1 and me.Zn != 0:
      rv = 1. - me.n
      for k in range(2, me.n + 1):
        rv *= 2./k
      rv *= me.Zn
      me.nr = 0
    else:
      me.nr = me.nrsm[1] / me.nrsm[0]
      rv = -me.Bring / me.nr;
    me.fb0 = me.fbsm[0][1] / me.fbsm[0][0]
    if me.tag == 'M': # original mcrat.c
      me.nz = me.nzsm[1] / me.nzsm[0]
      me.ta0 = me.tacc[0][1] / me.tacc[0][0]
      me.ta1 = me.tacc[1][1] / me.tacc[1][0]
      me.fb1 = (me.fbsm[0][1] + me.fbsm[1][1]) / (
                me.fbsm[0][0] + me.fbsm[1][0] / me.nz)
      me.fb2 = me.fbsm[2][1] / me.fbsm[2][0] * (me.Z[2]
               / me.Z[1]) * me.ta0 / me.ta1 * me.nz
      me.vir = [me.fb0 * rv, me.fb1 * rv, me.fb2 * rv]
    else: # mcrat0.c
      if me.use3:
        # vir, fb, nr
        me.vir = [me.fb0 * rv, me.fb0, me.nr]
      else:
        me.vir = [me.fb0 * rv]



  def clear(me):
    ''' me = 0 '''
    me.fbsm = [[0, 0],] * 3
    me.nzsm = [0, 0]
    me.nrsm = [0, 0]
    me.tacc = [[0, 0],] * 2
    me.tot = 0



  def absorb(a, b):
    ''' a += b '''
    for i in range(3):
      xinc(a.fbsm[i], b.fbsm[i])
    xinc(a.nzsm, b.nzsm)
    xinc(a.nrsm, b.nrsm)
    for i in range(2):
      xinc(a.tacc, b.tacc)
    a.tot += b.tot



  def aggregate(mr, ls):
    ''' mr = sum(ls) for int* '''
    mr.clear()
    for mr1 in ls:
      mr.absorb(mr1) # mr += mr1

    mr.recompute() # on all data

    # compute errors
    m = len(ls)
    virls = [None] * m
    for i in range(m): # loop over copies
      virls[i] = ls[i].vir
    mr.evir = errarr( virls )



  def save(me, fn):
    s = ""
    if me.tag == 'M': # output of mcrat (obsolete)
      info = "#M %d %d 1 %.14e %.14e V0\n" % (
          me.D, me.n, me.Z[1], me.Z[2])
      for k in range(3):
        s += "%16.0f %+17.14f " % (
            me.fbsm[k][0], me.fbsm[k][1] / me.fbsm[k][0])
      s += "%16.0f %+17.14f %16.0f %+18.14f " % (
            me.nzsm[0], me.nzsm[1] / me.nzsm[0],
            me.nrsm[0], me.nrsm[1] / me.nrsm[0])
      s += "%16.0f %16.14f %16.0f %16.14f " % (
            me.tacc[0][0], me.tacc[0][1] / me.tacc[0][0],
            me.tacc[1][0], me.tacc[1][1] / me.tacc[1][0])
      s += "%+20.14e %+20.14e %+20.14e " % (
            me.vir[0], me.vir[1], me.vir[2])
      if me.evir:
        s += "%9.2e %9.2e %9.2e " % (
            me.evir[0], me.evir[1], me.evir[2])
      if verbose:
        print "saved mr file %s, tot %g, %g, %g" % (
            fn, me.fbsm[0][0], me.fbsm[1][0], me.fbsm[2][0])
        print "D %d, n %d, %+.6e (%.1e) %+.6e (%.1e) %+.6e (%.1e) " % (me.D, me.n,
            me.vir[0], me.evir[0], me.vir[1], me.evir[1], me.vir[2], me.evir[2])
    else: # output of mcrat0
      if me.version >= 1:
        info = "#0 %d %d V1 %20.14e %20.14e\n" % (me.D, me.n, me.Bring, me.Zn)
      else:
        info = "#0 %d %d V0\n" % (me.D, me.n)
      s += "%16.0f %+17.14f " % (
            me.fbsm[0][0], me.fbsm[0][1] / me.fbsm[0][0])
      s += "%16.0f %+18.14f " % (
            me.nrsm[0], me.nrsm[1] / me.nrsm[0])
      s += "%+20.14e " % (me.vir[0])
      if me.evir:
        s += "%9.2e " % (me.evir[0])
      if verbose:
        print "D %d, n %d, %+.6e (%.1e), " % (me.D, me.n,
            me.vir[0], me.evir[0]),
        print "saved mr file %s, tot %g (%7.1e/%7.1e)" % (
            fn, me.tot, me.have, me.hsig)
    s += "\n"
    open(fn, "w").write(info + s)



  @staticmethod
  def sumdat(fnbas, fnls, use3 = 0):
    ''' sum over data in fnbas + fnls for mr*.dat '''
    mr = MR(fnbas, use3 = use3) # create a class, read data in fnbas
    mr.clear()
    m = len(fnls) + 1
    fnls = [fnbas,] + fnls
    mrls = []
    # sum over all file lists
    for i in range(m):
      fn = fnls[i]
      mr1 = MR(fn, use3 = use3) # read data in fn
      mr1.recompute() # compute individual vir
      mrls += [mr1, ]

    # histogram analysis
    totls = [li.tot for li in mrls]
    mr.have, mr.hsig, hisbad = varanalysis(
        totls, [li.fn for li in mrls],
        "histogram", 5.0, 0.1, verbose)

    excl = []
    while 1:
      fns = [li.fn for li in mrls]
      # variance analysis
      mr.virave, mr.virsig, virbad = varanalysis(
          [li.vir[0] for li in mrls], fns,
          "virial", mulsig, 0.0, verbose)

      if not rmbad: break
      badls = virbad
      if len(badls) >= len(mrls): # exclude everything
        print "the list is so bad that I have to exclude everything!"
        break
      excl += [fns[i] for i in badls]
      if len(badls) == 0: break
      # exclude the bad items
      m = len(mrls)
      mrls = [mrls[i] for i in range(m) if not i in badls]

    # aggregate after the exclusion
    if len(excl):
      print "excluded %s" % excl
    mr.aggregate(mrls)

    fnout = fnbas + "a"
    mr.save(fnout)
    return (mr.vir, mr.evir, mr.tot)





class GC:
  ''' a class for combining data from mcgc2.c and mcgcr2.c '''
  once = 0

  def __init__(me, fn, fnBring = None):
    s = open(fn).readlines()
    if s[0][1] in ('R', 'S'):
      me.initZrr(s[0])
      me.loadZrr(s[1:])
    elif s[0][1] == 'H':
      me.initZrh(s[0])
      me.loadZrh(s[1:])
    me.loadBring(fnBring)



  def loadBring(me, fn = None):
    ''' load virial coefficients from Bring.dat '''
    if not fn: fn = findBring()

    ls = open(fn).readlines()
    for s in ls:
      bs = s.split()
      if int(bs[0]) == me.D:
        break
    else:
      print "%s: no data for D %d" % (fn, me.D)
      return
    nmax = min(len(bs) - 1, me.nmax)
    for n in range(1, nmax + 1):
      me.Bring[n] = fabs(float(bs[n]))
    if not GC.once:
      if verbose:
        print "loaded %d values from %s" % (nmax, fn)
      GC.once = 1


  def clear(a):
    ''' a = 0 '''
    for i in range(a.nens):
      # accumulate data
      a.hist[i] = 0
      a.nup[i] = [1e-20, 0]
      a.ndown[i] = [1e-20, 0]
      if a.tag == 'H':
        a.ngpr[i] = [1e-20, 0]
        n = i
      else:
        n = a.n[i]
      if a.tag == 'H' or a.tp[i] == 0:
        a.ring[n] = [1e-20, 0]
        a.nedg[n] = [1e-20, 0]
        a.ncsp[n] = [1e-20, 0]
        a.fbsm[n] = [1e-20, 0, 0]
    a.tot = 0



  def absorb(a, b):
    ''' a += b '''
    for i in range(a.nens):
      # accumulate data
      a.hist[i] += b.hist[i]
      xinc(a.nup[i], b.nup[i])
      xinc(a.ndown[i], b.ndown[i])
      if a.tag == 'H':
        xinc(a.ngpr[i], b.ngpr[i])
        n = i
      else:
        n = a.n[i]
      if a.tag == 'H' or a.tp[i] == 0:
        xinc(a.ring[n], b.ring[n])
        xinc(a.nedg[n], b.nedg[n])
        xinc(a.ncsp[n], b.ncsp[n])
        xinc(a.fbsm[n], b.fbsm[n])
    a.tot += b.tot



  def recompute(me):
    if me.tag in ('R', 'S'):
      (me.Zr1, me.rc1, me.sr1) = me.updateZrr()
      me.computeZrr(me.Zr1, me.rc1, me.sr1)
    elif me.tag == 'H' or not me.tag:
      (me.Zr1, me.rc1) = me.updateZrh()
      me.computeZrh(me.Zr1, me.rc1)
    else:
      print "do not support %s" % me.tag
      raise Exception



  def initZrh(me, s):
    ''' initialize Zrh data structure '''
    me.info = s.split()
    me.tag = me.info[0][1:] # the letter after '#', if any
    me.D = int(me.info[1])
    me.nmax = int(me.info[2])
    me.nens = me.nmax + 1
    me.sver = me.info[3]
    me.ver = getver(me.sver)
    me.nmin = int(me.info[4])
    me.nedxmax = me.info[-1]
    cnt = me.nmax + 1
    me.Zr = [0] * cnt
    me.rc = [0] * cnt
    me.Z = [1] * cnt
    me.hist = [0] * cnt
    me.tot = 0
    me.nup = [[1e-20, 0],] * cnt
    me.ndown = [[1e-20, 0],] * cnt
    me.ngpr = [[1e-20, 0],] * cnt
    me.ring = [[1e-20, 0],] * cnt
    me.nedg = [[1e-20, 0],] * cnt
    me.ncsp = [[1e-20, 0],] * cnt
    me.fbsm = [[1e-20, 0, 0],] * cnt
    me.B = [0] * cnt
    me.B2 = [0] * cnt
    me.Bring = [1] * cnt
    me.Zn = [1] * cnt
    me.ZZ = [1] * cnt
    me.eB = me.eB2 = me.eZn = me.eZZ = None



  def loadZrh(me, s):
    ''' absorb data in me '''
    me.tot = 0
    for ln in s:
      # parse the lines
      x = ln.strip().split()
      i = int(x[0])
      # parameters
      me.ZZ[i] = float(x[1])
      me.Z[i] = float(x[2])
      me.Zr[i] = float(x[3])
      me.rc[i] = float(x[4])
      # accumulate data
      me.hist[i] = float(x[5])
      me.tot += me.hist[i]
      me.nup[i] = toarr((x[6], x[8]))
      me.ndown[i] = toarr((x[7], x[9]))
      me.ngpr[i] = toarr((x[10], x[11]))
      k = 12
      if me.ver >= 4: # load ring data
        me.ring[i] = toarr((x[k], x[k+1]))
        k += 2
      me.nedg[i] = toarr((x[k], x[k+3]))
      me.ncsp[i] = toarr((x[k+1], x[k+4]))
      me.fbsm[i] = toarr((x[k+2], x[k+6], x[k+5]))



  def updateZrh(me, mindata = 100):
    ''' update Zrh file '''
    Zr = me.Zr[:]
    rc = me.rc[:]
    for i in range(1, me.nens):
      r = 1
      if me.nup[i][1] > mindata and me.ndown[i][1] > mindata:
        r = me.nup[i][1] / me.nup[i][0] / (me.ndown[i][1] / me.ndown[i][0])
      if me.tag == 'H' and me.ngpr[i][1] > mindata:
        nZr = me.ngpr[i][0] / me.ngpr[i][1]
        rc[i] *= pow(r * Zr[i] / nZr, 1./me.D)
        Zr[i] = nZr
      else:
        Zr[i] *= r
      if not me.tag:
        rc[i] = pow(Zr[i], 1./me.D)
    return (Zr, rc)



  def computeZrh(me, Zr, rc):
    ''' compute the virial series '''
    z = 1
    fac = 1
    prevZ = 1
    for n in range(1, me.nmax + 1):
      z *= Zr[n]
      if me.tag == 'H':
        z *= pow(rc[n], me.D) * n * (n - 1)
      if n <= 2:
        z = 1
      elif n == 3 or (n == 4 and me.D <= 12):
        z = me.Z[n]
      me.Zn[n] = me.Z[n] = z
      if n >= 2: fac *= 2./n
      if n <= 2:
        fbav = 1
      else:
        fbav = me.fbsm[n][1] / me.fbsm[n][0]
      nr = me.ring[n][1] / me.ring[n][0]
      me.B[n] = (1. - n) * fac * z * fbav
      if nr > 0:
        me.B2[n] = -fbav / nr * me.Bring[n]
      me.ZZ[n] = me.Zn[n] / prevZ
      prevZ = me.Zn[n]



  def saveZrh(me, fn, alt = False):
    ''' return a string of the compact Zrh file '''
    s = [""] * (me.nmax + 1)
    (Zr, rc) = (me.Zr, me.rc)
    if alt:
      (Zr, rc) = (me.Zr1, me.rc1)
    for n in range(1, me.nens):
      s[n] = "%3d %20.14e %20.14e %18.14f %18.14f " % (
          n, me.ZZ[n], me.Zn[n], Zr[n], rc[n])
      s[n] += "%14.0f %14.0f %14.0f %.14f %.14f " % (
          me.hist[n], me.nup[n][0], me.ndown[n][0],
          me.nup[n][1] / me.nup[n][0],
          me.ndown[n][1] / me.ndown[n][0])
      if me.tag == 'H':
        s[n] += "%14.0f %17.14f " % (
          me.ngpr[n][0], me.ngpr[n][1] / me.ngpr[n][0])
      if me.ver >= 4:
        s[n] += "%14.0f %20.14e " % (
            me.ring[n][0], me.ring[n][1] / me.ring[n][0])
      s[n] += "%14.0f %14.0f %14.0f " % (
          me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
      s[n] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
          me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
          me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
          me.B[n])
      if me.ver >= 4:
        s[n] += "%+20.14e " % me.B2[n]
      if me.eB:
        if me.ver >= 4:
          s[n] += "%9.2e " % me.eB2[n]
        s[n] += "%9.2e %9.2e %9.2e " % (
            me.eB[n], me.eZZ[n], me.eZn[n])
      s[n] = s[n].rstrip() + "\n"
    s[0] = "#%s %d %d V%d %d 1 %s\n" % (me.tag, me.D,
        me.nmax, me.ver, me.nmin, me.nedxmax)
    open(fn, "w").writelines(s)
    if verbose:
      print "saved Zrh file %s, tot %g" % (fn, me.tot)



  def initZrr(me, s):
    ''' initialize Zrr data structure '''
    me.info = s.split()
    me.tag = me.info[0][1:] # the letter after '#', if any
    me.D = int(me.info[1])
    me.nens = int(me.info[2])
    me.ens0 = int(me.info[3])
    me.nmin = int(me.info[4])
    me.nmax = int(me.info[5])
    me.m = int(me.info[6])
    me.sver = me.info[7]
    me.ver = getver(me.sver)
    me.nedxmax = me.info[-1]
    me.tp = (range(me.m) * (me.nmax + 1))[:me.nens]
    me.n = [ int((q + me.m - 0.999)/me.m) for q in range(me.nens) ]
    me.Zr = [0] * me.nens
    me.rc = [0] * me.nens
    me.sr = [0] * me.nens
    me.Z = [1] * me.nens
    me.hist = [0] * me.nens
    me.tot = 0
    me.nup = [[1e-20, 0],] * me.nens
    me.ndown = [[1e-20, 0],] * me.nens
    cnt = me.nmax + 1
    me.ring = [[1e-20, 0],] * cnt
    me.nedg = [[1e-20, 0],] * cnt
    me.ncsp = [[1e-20, 0],] * cnt
    me.fbsm = [[1e-20, 0, 0],] * cnt
    me.B = [0] * cnt
    me.B2 = [0] * cnt
    me.Bring = [0] * cnt
    me.Zn = [1] * cnt
    me.ZZ = [1] * cnt



  def loadZrr(me, s):
    ''' load data to me '''
    me.tot = 0
    for ln in s:
      # parse the lines
      x = ln.strip().split()
      i = int(x[0])
      # parameters
      me.tp[i] = int(x[1])
      me.n[i] = int(x[2])
      me.Zr[i] = float(x[3])
      me.Z[i] = float(x[4])
      me.rc[i] = float(x[5])
      me.sr[i] = float(x[6])
      # accumulate data
      me.hist[i] = float(x[7])
      me.tot += me.hist[i]
      me.nup[i] = toarr((x[8], x[10]))
      me.ndown[i] = toarr((x[9], x[11]))
      if me.tp[i] == 0:
        n = me.n[i]
        k = 12
        if me.ver >= 4:
          me.ring[n] = toarr((x[k], x[k+1]))
          k += 2
        me.nedg[n] = toarr((x[k], x[k+3]))
        me.ncsp[n] = toarr((x[k+1], x[k+4]))
        me.fbsm[n] = toarr((x[k+2], x[k+6], x[k+5]))



  def updateZrr(me, mindata = 100):
    Zr = me.Zr[:]
    rc = me.rc[:]
    sr = me.sr[:]
    for i in range(me.ens0 + 1, me.nens):
      if me.nup[i][1] > mindata and me.ndown[i][1] > mindata:
        r = me.nup[i][1] / me.nup[i][0] / (me.ndown[i][1] / me.ndown[i][0])
      else:
        r = 1
      if me.tp[i - 1] == 0:
        Zr[i] *= r
      else:
        sr[i - 1] *= pow(r, 1./me.D)
    return (Zr, rc, sr)



  def computeZrr(me, Zr, rc, sr):
    ''' compute the virial series '''
    z = 1
    fac = 1
    prevZ = 1
    for i in range(1, me.nens):
      n = me.n[i]
      z *= Zr[i]
      if me.tp[i - 1] != 0:
        z *= pow(sr[i - 1], me.D)
      else:
        z *= pow(rc[i], me.D)
      if me.tp[i] == 0:
        if me.tag == 'R':
          z *= n * (n - 1) * .5
        if n <= 2:
          z = 1
        elif n == 3 or (n == 4 and me.D <= 12):
          z = me.Z[i]
      me.Z[i] = z
      if me.tp[i] == 0:
        if n >= 2: fac *= 2./n
        if n <= 2:
          fbav = 1
        else:
          fbav = me.fbsm[n][1] / me.fbsm[n][0]
        nr = me.ring[n][1] / me.ring[n][0]
        me.B[n] = (1. - n) * fac * z * fbav
        if nr > 0:
          me.B2[n] = -fbav / nr * me.Bring[n]
        me.Zn[n] = me.Z[i]
        me.ZZ[n] = me.Zn[n] / prevZ
        prevZ = me.Zn[n]



  def saveZrr(me, fn, alt = False):
    ''' output to Zrr file '''
    s = [""] * (me.nens + 1)
    (Zr, rc, sr) = (me.Zr, me.rc, me.sr)
    if alt: # alternative buffer
      (Zr, rc, sr) = (me.Zr1, me.rc1, me.sr1)
    for i in range(1, me.nens):
      n = me.n[i]
      tp = me.tp[i]
      s[i] = "%4d %d %3d %20.14e %20.14e %18.14f %18.14f " % (
          i, me.tp[i], me.n[i], Zr[i], me.Z[i], rc[i], sr[i])
      s[i] += "%14.0f %14.0f %14.0f %.14f %.14f " % (
          me.hist[i], me.nup[i][0], me.ndown[i][0],
          me.nup[i][1] / me.nup[i][0],
          me.ndown[i][1] / me.ndown[i][0])
      if tp == 0:
        if me.ver >= 4:
          s[i] += "%14.0f %17.14f " % (
              me.ring[n][0], me.ring[n][1] / me.ring[n][0])
        s[i] += "%14.0f %14.0f %14.0f " % (
            me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
        s[i] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
            me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
            me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
            me.B[n])
        if me.ver >= 4:
          s[i] += "%+20.14e " % me.B2[n]
        if me.eB:
          if me.ver >= 4:
            s[i] += "%9.2e " % me.eB2[n]
          s[i] += "%9.2e %9.2e %9.2e " % (
              me.eB[n], me.eZZ[n], me.eZn[n])
      s[i] = s[i].rstrip() + "\n"
    s[0] = ' '.join(me.info) + '\n'
    open(fn, "w").writelines(s)
    if verbose:
      print "saved Zrr file %s, tot %g" % (fn, me.tot)



  def saveZr(me, fn):
    ''' return a string of the compact Zr file '''
    s = [""] * (me.nmax + 1)
    for i in range(1, me.nens):
      tp = me.tp[i]
      if tp != 0: continue
      n = me.n[i]
      s[n] = "%3d %18.14f %20.14e " % (
          n, me.ZZ[n], me.Zn[n])
      s[n] += "%14.0f %14.0f %14.0f %.14f %.14f " % (
          me.hist[i], me.nup[i][0], me.ndown[i][0],
          me.nup[i][1] / me.nup[i][0],
          me.ndown[i][1] / me.ndown[i][0])
      if me.ver >= 4:
        s[n] += "%14.0f %17.14f " % (
            me.ring[n][0], me.ring[n][1] / me.ring[n][0])
      s[n] += "%14.0f %14.0f %14.0f " % (
          me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
      s[n] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
          me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
          me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
          me.B[n])
      if me.ver >= 4:
        s[n] += "%+20.14e " % me.B2[n]
      if me.eB:
        if me.ver >= 4:
          s[n] += "%9.2e " % me.eB2[n]
        s[n] += "%9.2e %9.2e %9.2e " % (
            me.eB[n], me.eZZ[n], me.eZn[n])
      s[n] = s[n].rstrip() + "\n"
    s[0] = "# %d %d V3 %d 1 %s\n" % (
        me.D, me.nmax, me.nmin, me.nedxmax)
    open(fn, "w").writelines(s)
    if verbose:
      print "saved Zr file %s, tot %g" % (fn, me.tot)



  @staticmethod
  def sumdat(fnbas, fnls, fnout = None, fnZr = None):
    ''' sum over data in fnbas + fnls, output to fnZrr '''
    gc = GC(fnbas)
    m = len(fnls) + 1
    gcls = [None] * m
    fnls = [fnbas,] + fnls
    # sum over all file lists
    for i in range(m):
      fn = fnls[i]
      gcls[i] = GC(fn)
      gcls[i].recompute()
      if i > 0: gc.absorb(gcls[i])
    gc.recompute()

    # compute errors
    Bls  = [None] * m
    B2ls = [None] * m
    Znls = [None] * m
    ZZls = [None] * m
    for i in range(m):
      Bls[i]  = gcls[i].B
      B2ls[i] = gcls[i].B2
      Znls[i] = gcls[i].Zn
      ZZls[i] = gcls[i].ZZ
    gc.eB  = errarr( Bls )
    gc.eB2 = errarr( B2ls )
    gc.eZn = errarr( Znls )
    gc.eZZ = errarr( ZZls )

    if not fnout: fnout = fnbas + "a"
    if gc.tag in ('R', 'S'):
      gc.saveZrr(fnout)
      # refined parameters
      gc.saveZrr("refine.data", True)
      fnZr = "Zr" + fnout
      if fnout.startswith("Zrr"):
        fnZr = "Zr" + fnout[3:]
      gc.saveZr(fnZr)
    else:
      gc.saveZrh(fnout)
      gc.saveZrh("refine.data", True)
    return gc.B2, gc.eB2, gc.tot




def niceprint(x, err, n = 0, i = 0, strtot = ""):
  ''' print the number with error aligned on the second line '''
  ox = floor(log10(fabs(x)))
  oy = floor(log10(err))
  if ox >= oy:
    padding = int(ox - oy + .5)
  else:
    padding = int(ox - oy - .5)
  if padding <= 0: padding -= 1
  if verbose:
    print "%4d/%4d: %+20.10e   %s\n" % (n, i, x, strtot) + " " * (
        padding + 15) + "%9.2e" % (err)



def getnstfb(d):
  ''' guess the interval of computing the star and ring contents
      from the name of the directory d '''
  m = re.search("w([0-9]*)", d)
  if m == None: return 1
  return int( m.group(1) )



def getlserr(ls, n, usetot = 0):
  ''' combine data from different simulations
      `usetot': use total instead of the inverse variance
                as the relative weight
      return average, error, total, strtot
  '''
  nsim = len(ls) # number of simulations
  if nsim <= 0: return None, None, None, None

  nq = len(ls[0][0]) # number of quantities
  if verbose:
    print "list has %s simulations, each with %d quantities" % (nsim, nq)

  # construct a dictionary of totals of different nstfb
  totdic = {}
  tot = 0
  for sim in range(nsim): # loop over simulations
    # guess nstfb
    nstfb = getnstfb(ls[sim][3])
    # update the total
    if nstfb in totdic:
      totdic[nstfb] += ls[sim][2]
    else:
      totdic[nstfb] = ls[sim][2]
    tot += ls[sim][2]

  # construct a string of totals
  strtot = "total: "
  for nstfb in totdic:
    strtot += "%.1e(%d), " % (totdic[nstfb], nstfb)
  strtot = strtot.strip()

  avls = []
  errls = []
  for q in range(nq): # loop over quantities, usually just one quanity
    suminvvar = 0
    sumx = 0
    sumw = 1e-100
    sumvar = 0
    # ls[sim][0][q] is the average
    # ls[sim][1][q] is the standard deviation
    # ls[sim][2] is the total
    wt = [0] * nsim
    for sim in range(nsim): # loop over simulations
      if ls[sim][0][q] == 0: # missing e.g. the ring content
        continue
      if usetot:
        wt[sim] = ls[sim][2]
        sumvar += (ls[sim][1][q] * wt[sim])**2
      else:
        invsig = 1./ls[sim][1][q]
        invvar = invsig * invsig # inverse variance
        wt[sim] = invvar
        suminvvar += invvar
      sumw += wt[sim]
      sumx += ls[sim][0][q] * wt[sim]
    if nsim > 1 and verbose: # print the break down
      for sim in range(nsim):
        if ls[sim][0][q] == 0: # missing e.g. the ring content
          continue
        # id, dir, weight, x, err
        print "%2d %-30s %5.2f%% %+20.10e %9.2e %16.0f" % (
            sim + 1, ls[sim][3], 100.*wt[sim]/sumw,
            ls[sim][0][q], ls[sim][1][q], ls[sim][2])
    avx = sumx / sumw
    if usetot:
      err = sqrt(sumvar) / sumw
    else:
      err = 1./sqrt(suminvvar)
    if avx == 0 and err == 0:
      #print "skip: ave %s, err %s, n %d" % (avx, err, n)
      pass
    else:
      niceprint(avx, err, n, q, strtot)
    avls += [avx,]
    errls += [err,]
  return avls, errls, tot, strtot



def guessdirs(n):
  # try to search directories of order n
  dirs = [d for d in glob.glob("n%d*" % n) if os.path.isdir(d)]
  dirs += [d for d in glob.glob("n%d*/mic*" % n) if os.path.isdir(d)]
  dirs += [d for d in glob.glob("n%d*/n%dmic*" % (n, n)) if os.path.isdir(d)]
  dirs = list( set( dirs ) )
  dirs.sort()
  if len(dirs) == 0:
    print "cannot find directories of order %d" % n
    return None
  return dirs



def dodirs(dirs, n, sum3 = 0):
  ''' handle data for the same order '''
  if dirs == None or len(dirs) == 0:
    dirs = guessdirs(n)
    if dirs == None:
      return None, None, None, None
  # dictionary of list with different nstfb,
  # the interval of computing the star and ring contents
  dicls = {}
  i = 0
  for d in dirs:
    # get nstfb from the directory name
    nstfb = getnstfb(d)
    if nstfb in dicls: ls = dicls[nstfb]
    else: ls = []
    (x, err, tot) = sumdat_wrapper(None, d, sum3 = sum3)
    if x:
      # average (array), standard deviation(array), total, directory, name
      ls += [(x, err, tot, d),]
      i += 1
      if verbose:
        print "Directory %d: %s, %s (%s), tot %s\n" % (i, d, x, err, tot)
    dicls[nstfb] = ls

  ''' 1. we first aggregate directories with the same nstfb
         with the total being the weights
      2. we then combine data with different nstfb '''
  newls = []
  strtot = ""
  for nstfb in dicls:
    ls = dicls[nstfb]
    # combine data directories with the same nstfb
    av, err, tot, strtot = getlserr(ls, n, usetot = 1)
    if av == None: continue;
    newls += [ (av, err, tot, "w" + str(nstfb)), ]
  # combine data with different nstfb
  return getlserr(newls, n, usetot = 0)



def scanorders(fns):
  ''' scan over orders '''
  n = 4
  ls = []
  while 1:
    ret = dodirs(None, n)
    #print ret
    #raw_input()
    x, err, tot, strtot = ret[0], ret[1], ret[2], ret[3]
    if x == None: break
    ls += [(n, x, err, strtot)]
    n += 1
  for n, x, err, strtot in ls:
    for i in range(len(x)):
      niceprint(x[i], err[i], n, i, strtot)
  return ls



def tofnbas(fn):
  return os.path.splitext(fn)[0] + ".dat"



def getfnbas(fninp):
  ''' get the input file '''
  if fninp:
    return tofnbas(fninp)

  # guess the input file name
  # we only need Zrr files, but not Zr files
  ls = glob.glob("ZrrD*n*.dat[1-9]*")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("ZrhD*n*.dat[1-9]*")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("mrD*n*.dat[1-9]*")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("intD*n*.dat[1-9]*")
  if len(ls): return tofnbas(ls[0])
  return None



def getslaves(fnbas):
  ''' get a list file names of non-master nodes '''
  ls = glob.glob(fnbas + "[0-9]*")
  if verbose:
    print "got %d files for %s" % (len(ls) + 1, fnbas)
  return ls



def checkdirfn(dir, fn):
  ''' check if the directory name and file name fn represent
      the same simulation '''
  # check dimensions
  mdir = re.search(r"D([0-9]+)", dir)
  if mdir:
    dird = int( mdir.group(1) )
    mfn = re.search(r"D([0-9]+)", fn)
    if mfn:
      fnd = int( mfn.group(1) )
      if dirnamecheck and dird != fnd:
        print "file name %s (D = %d) incompatible with the dir name %s (D = %d)" % (
            fn, fnd, dir, dird)
        return 1

  # check orders
  mdir = re.search("n([0-9]+)mic", dir)
  if not mdir:
    mdir = re.search("n([0-9]+)", dir)
  if mdir:
    dirn = int( mdir.group(1) )
    mfn = re.search("n([0-9]+)", fn)
    if mfn:
      fnn = int( mfn.group(1) )
      if dirnamecheck and dirn != fnn:
        print "file name %s (n = %d) incompatible with the dir name %s (n = %d)" % (
            fn, fnn, dir, dirn)
        return 1

  return 0



def sumdat_wrapper(fninp, dir = None, sum3 = 0, checkdirname = True):
  ''' summarize data for fninp under a single directory '''

  dir0 = os.getcwd()
  if dir: os.chdir(dir)
  # guess the basic file name if necessary
  fnbas = getfnbas(fninp)
  if fnbas:
    # check if the file name is compatible with directory name
    if checkdirname and checkdirfn(os.getcwd(), fnbas) != 0:
      if dir:
        os.chdir(dir0)
      return None, None, None
  if not fnbas: # no basic file name
    if verbose:
      print "cannot find a .dat file under %s" % dir
    if dir:
      os.chdir(dir0)
    return None, None, None
  elif not os.path.exists(fnbas):
    if verbose:
      print "cannot find %s, dir %s %s" % (fnbas, dir, dir0)
    if dir:
      os.chdir(dir0)
    return None, None, None

  # obtain slave files from the basic file
  fnls = getslaves(fnbas)

  if fnbas.startswith("Z"):
    tp = "Z"
  elif fnbas.startswith("mr"):
    tp = "mr"
  elif fnbas.startswith("intg"):
    tp = "intg"
  else:
    sinfo = open(fnbas).readlines()[0]
    if sinfo[1] in "RS":
      tp = "Z"
    else:
      tp = "mr"

  if tp == "Z":
    (x, err, tot) = GC.sumdat(fnbas, fnls)
  elif tp == "mr":
    (x, err, tot) = MR.sumdat(fnbas, fnls, sum3)
  else:
    (x, err, tot) = INT.sumdat(fnbas, fnls)
  if dir: os.chdir(dir0)
  return x, err, tot



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] [input.dat]"
  print """
  Aggregate multiple MPI/OpenMP data generated by
    hsrh.c, mcrat0.c, mcgc2.c, mcgcr2.c
  The input are something like
    intD3n6.dat, mrD3n9.dat, ZrhD8n64.dat, ZrrD100r4n64.dat
  Without the input, the program will try search proper input dat files

  OPTIONS:
   -s:     scan over directories of different n
   -x:     exclude data directories with large deviations from the average
   -y[3]:  multiples of the standard deviation for the above exclusion
   -n:     grand aggregation for all folders with n
   -d:     input are directories to aggregate
   -t:     use the total instead of the inverse variance to weight individual
           simulations during the above aggregation
   -N:     no check for the directory name

  EXAMPLES:
  1. aggregrate data of mrD2n10.dat, mrD2n10.dat1, ...
    $PROG  mrD2n10.dat

  2. aggregate directories dir1, dir2, ... of order n = 10
    $PROG -d -n 10 dir1 dir2 dir3

  3. search directories for order 10 and aggregate them
    $PROG -n 10

  4. scan directories of all orders and gregrates them by orders
    $PROG -s

  """.replace("$PROG", sys.argv[0])
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "sdn:xy:N",
        ["scan", "dir", "rmbad", "order=", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, rmbad, mulsig, dirnamecheck
  n = 0
  isdir = 0

  for o, a in opts:
    if o in ("-d", "--dir"):
      isdir = 1
    elif o in ("-s", "--scan"):
      isdir = -1
    elif o in ("-n", "--order",):
      n = int(a)
      isdir = 1
    elif o in ("-x", "--rmbad"):
      rmbad = 1
    elif o in ("-y",):
      mulsig = float(a)
    elif o in ("-N",):
      dirnamecheck = False
    elif o in ("--verbose=",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  fns = args
  if isdir > 0:
    if n == 0: usage()
  return isdir, fns, n



# main function starts here
if __name__ == "__main__":
  isdir, fns, n = doargs()

  if isdir < 0: # scan different orders
    scanorders(fns)
  elif isdir > 0: # multiple directory mode
    dodirs(fns, n)
  else: # under a single directory
    if len(fns):
      for fn in fns: sumdat_wrapper(fn)
    else:
      sumdat_wrapper(None, checkdirname = False)

