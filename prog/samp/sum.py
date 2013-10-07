#!/usr/bin/env python

import os, sys, re, glob
from math import *



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



class INT:
  ''' a class for combining data from hsrh.c '''
  once = 0

  def __init__(me, fn):
    s = open(fn).readlines()
    me.initint(s[0])
    me.adddata(s[1:])



  def initint(me, s):
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


  def adddata(me, s):
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
    ''' sum over data in fnbas + fnls '''
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
    intg.recompute()

    # compute errors
    virls = [None] * nls
    for i in range(nls): # loop over copies
      virls[i] = [intgls[i].vir,]
    intg.evir = errarr( virls )[0]

    fnout = fnbas + "a"
    intg.save(fnout)
    print "saving %s, vir %.8e, %+9.2e" % (fnout, intg.vir, intg.evir)



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



def findBring():
  ''' find Bring.dat '''
  # current directory
  fn = fn0 = "Bring.dat"
  if os.path.exists(fn): return fn
  # parent directory
  fn = os.path.join(os.pardir, fn0)
  if os.path.exists(fn): return fn
  # directory of the script
  d = os.path.split(os.path.realpath(__file__))[0]
  fn = os.path.join(d, fn0)
  if os.path.exists(fn): return fn
  print "Error: cannot load Bring data"
  return None



class MR:
  ''' a class for combining data from mcrat.c '''
  once = 0

  def __init__(me, fn, fnBring = None):
    s = open(fn).readlines()
    me.initmr(s[0])
    me.adddata(s[1])
    me.loadBring(fnBring)



  def loadBring(me, fn):
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
    me.Bring = fabs( float(bs[ me.n ]) )
    if not GC.once:
      print "loaded n %d, D %d from %s, %s" % (me.n, me.D, fn, me.Bring)
      GC.once = 1



  def initmr(me, s):
    arr = s.split()
    me.tag = arr[0][1:].strip()
    me.D = int(arr[1])
    me.n = int(arr[2])
    if me.tag == 'M':
      me.Z = [1, float(arr[4]), float(arr[5])]
    me.fbsm = [[0,0], [0,0], [0,0]]
    me.nzsm = [0,0]
    me.nrsm = [0,0]
    me.tacc = [[0,0], [0,0]]
    me.evir = None



  def adddata(me, s):
    x = s.strip().split()
    me.fbsm[0] = toarr((x[0], x[1]))
    if me.tag == 'M':
      me.fbsm[1] = toarr((x[2], x[3]))
      me.fbsm[2] = toarr((x[4], x[5]))
      me.nzsm = toarr((x[6], x[7]))
      me.nrsm = toarr((x[8], x[9]))
      me.tacc[0] = toarr((x[10], x[11]))
      me.tacc[1] = toarr((x[12], x[13]))
    else:
      me.nrsm = toarr((x[2], x[3]))



  def absorb(a, b):
    ''' a += b '''
    for i in range(3):
      xinc(a.fbsm[i], b.fbsm[i])
    xinc(a.nzsm, b.nzsm)
    xinc(a.nrsm, b.nrsm)
    for i in range(2):
      xinc(a.tacc, b.tacc)



  def recompute(me):
    ''' recompute Zr '''
    me.nr = me.nrsm[1] / me.nrsm[0]
    rv = -me.Bring / me.nr;
    me.fb0 = me.fbsm[0][1] / me.fbsm[0][0]
    if me.tag == 'M':
      me.nz = me.nzsm[1] / me.nzsm[0]
      me.ta0 = me.tacc[0][1] / me.tacc[0][0]
      me.ta1 = me.tacc[1][1] / me.tacc[1][0]
      me.fb1 = (me.fbsm[0][1] + me.fbsm[1][1]) / (
                me.fbsm[0][0] + me.fbsm[1][0] / me.nz)
      me.fb2 = me.fbsm[2][1] / me.fbsm[2][0] * (me.Z[2]
               / me.Z[1]) * me.ta0 / me.ta1 * me.nz
      me.vir = [me.fb0 * rv, me.fb1 * rv, me.fb2 * rv]
    else:
      me.vir = [me.fb0 * rv]



  def save(me, fn):
    s = ""
    if me.tag == 'M': # output of mcrat
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
      print "saved mr file %s, tot %g, %g, %g" % (
          fn, me.fbsm[0][0], me.fbsm[1][0], me.fbsm[2][0])
      print "D %d, n %d, %+.6e (%.1e) %+.6e (%.1e) %+.6e (%.1e) " % (me.D, me.n,
          me.vir[0], me.evir[0], me.vir[1], me.evir[1], me.vir[2], me.evir[2])
    else: # output of mcrat0
      info = "#0 %d %d V0\n" % (me.D, me.n)
      s += "%16.0f %+17.14f " % (
            me.fbsm[0][0], me.fbsm[0][1] / me.fbsm[0][0])
      s += "%16.0f %+18.14f " % (
            me.nrsm[0], me.nrsm[1] / me.nrsm[0])
      s += "%+20.14e " % (me.vir[0])
      if me.evir:
        s += "%9.2e " % (me.evir[0])
      print "saved mr file %s, tot %g" % (fn, me.fbsm[0][0])
      print "D %d, n %d, %+.6e (%.1e) " % (me.D, me.n,
          me.vir[0], me.evir[0])
    s += "\n"
    open(fn, "w").write(info + s)



  @staticmethod
  def sumdat(fnbas, fnls):
    ''' sum over data in fnbas + fnls '''
    mr = MR(fnbas)
    m = len(fnls) + 1
    mrls = [None] * m
    fnls = [fnbas,] + fnls
    # sum over all file lists
    for i in range(m):
      fn = fnls[i]
      mrls[i] = MR(fn)
      mrls[i].recompute() # compute individual vir
      if i > 0: mr.absorb(mrls[i])
    mr.recompute()

    # compute errors
    virls = [None] * m
    for i in range(m): # loop over copies
      virls[i] = mrls[i].vir
    mr.evir = errarr( virls )

    fnout = fnbas + "a"
    mr.save(fnout)





class GC:
  ''' a class for combining data from mcgc2.c and mcgcr2.c '''
  once = 0

  def __init__(me, fn, fnBring = None):
    s = open(fn).readlines()
    if s[0][1] in ('R', 'S'):
      me.initZrr(s[0])
      me.adddataZrr(s[1:])
    elif s[0][1] == 'H':
      me.initZrh(s[0])
      me.adddataZrh(s[1:])
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
      print "loaded %d values from %s" % (nmax, fn)
      GC.once = 1



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



  def adddataZrh(me, s):
    ''' absorb data in me '''
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
    me.tot = 0
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
      me.tot += me.hist[n]
    s[0] = "#%s %d %d V%d %d 1 %s\n" % (me.tag, me.D,
        me.nmax, me.ver, me.nmin, me.nedxmax)
    open(fn, "w").writelines(s)
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



  def adddataZrr(me, s):
    ''' absorb data in me '''
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
    me.tot = 0
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
      me.tot += me.hist[i]
    s[0] = ' '.join(me.info) + '\n'
    open(fn, "w").writelines(s)
    print "saved Zrr file %s, tot %g" % (fn, me.tot)



  def saveZr(me, fn):
    ''' return a string of the compact Zr file '''
    s = [""] * (me.nmax + 1)
    me.tot = 0
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
      me.tot += me.hist[i]
    s[0] = "# %d %d V3 %d 1 %s\n" % (
        me.D, me.nmax, me.nmin, me.nedxmax)
    open(fn, "w").writelines(s)
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



def tofnbas(fn):
  return os.path.splitext(fn)[0] + ".dat"



def getfnbas():
  ''' get the input file '''
  if len(sys.argv) > 1:
    return tofnbas(sys.argv[1])
  # we only need Zrr files, but not Zr files
  ls = glob.glob("ZrrD*n*.dat1")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("ZrhD*n*.dat1")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("mrD*n*.dat1")
  if len(ls): return tofnbas(ls[0])
  ls = glob.glob("intD*n*.dat1")
  if len(ls): return tofnbas(ls[0])
  return None



def getslaves(fnbas):
  ''' get a list file names of non-master nodes '''
  ls = []
  i = 1
  while 1:
    fn = "%s%d" % (fnbas, i)
    if not os.path.exists(fn):
      break
    ls += [fn,]
    i += 1
  return ls



# main function starts here
if __name__ == "__main__":
  fnbas = getfnbas()
  if not fnbas:
    print "please enter a .dat file"
    raise Exception
  elif not os.path.exists(fnbas):
    print "cannot find %s" % fnbas
    raise Exception
  fnls = getslaves(fnbas)
  if fnbas.startswith("Z"):
    GC.sumdat(fnbas, fnls)
  elif fnbas.startswith("mr"):
    MR.sumdat(fnbas, fnls)
  else:
    INT.sumdat(fnbas, fnls)

