#!/usr/bin/env python

import os, sys, re, glob
from math import *


def toarr(b):
  a = [float(x) for x in b]
  if a[0] <= 0: a[0] = 1e-20
  for i in range(1, len(a)):
    a[i] *= a[0]
  return a



def xinc(a, b):
  for i in range(len(b)):
    a[i] += b[i]



def errarr(arr):
  ''' calculate the error of the array `arr' '''
  m = len(arr) # number of copies
  if m <= 1:
    print arr
    raise Exception
  n = len(arr[0])
  err = [0] * n
  for i in range(n): # loop over array
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



class GC:
  def __init__(me, fn):
    me.type = "Zr"
    s = open(fn).readlines()
    if s[0][1] in ('R', 'S'):
      me.initZrr(s[0])
      me.adddataZrr(s[1:])
    elif s[0][1] == 'H':
      me.initZrh(s[0])
      me.adddataZrh(s[1:])
    me.haserr = False



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
    me.nmin = int(me.info[4])
    me.nedxmax = me.info[-1]
    cnt = me.nmax + 1
    me.Zr = [0] * cnt
    me.rc = [0] * cnt
    me.Z = [1] * cnt
    me.hist = [0] * cnt
    me.nup = [[0,0],] * cnt
    me.ndown = [[0,0],] * cnt
    me.ngpr = [[0,0],] * cnt
    me.nedg = [[0,0],] * cnt
    me.ncsp = [[0,0],] * cnt
    me.fbsm = [[0,0,0],] * cnt
    me.B = [0] * cnt
    me.Zn = [1] * cnt
    me.ZZ = [1] * cnt



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
      me.nedg[i] = toarr((x[12], x[15]))
      me.ncsp[i] = toarr((x[13], x[16]))
      me.fbsm[i] = toarr((x[14], x[18], x[17]))



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
      me.B[n] = (1. - n) * fac * z * fbav
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
        s[n] += "%14.0f %.14f " % (
          me.ngpr[n][0], me.ngpr[n][1] / me.ngpr[n][0])
      s[n] += "%14.0f %14.0f %14.0f " % (
          me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
      s[n] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
          me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
          me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
          me.B[n])
      if me.haserr:
        s[n] += "%9.2e %9.2e %9.2e " % (
            me.eB[n], me.eZZ[n], me.eZn[n])
      s[n] = s[n].rstrip() + "\n"
      me.tot += me.hist[n]
    s[0] = "#%s %d %d V3 %d 1 %s\n" % (
        me.tag, me.D, me.nmax, me.nmin, me.nedxmax)
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
    me.nedxmax = me.info[-1]
    me.tp = (range(me.m) * (me.nmax + 1))[:me.nens]
    me.n = [ int((q + me.m - 0.999)/me.m) for q in range(me.nens) ]
    me.Zr = [0] * me.nens
    me.rc = [0] * me.nens
    me.sr = [0] * me.nens
    me.Z = [1] * me.nens
    me.hist = [0] * me.nens
    me.nup = [[0,0],] * me.nens
    me.ndown = [[0,0],] * me.nens
    cnt = me.nmax + 1
    me.nedg = [[0,0],] * cnt
    me.ncsp = [[0,0],] * cnt
    me.fbsm = [[0,0,0],] * cnt
    me.B = [0] * cnt
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
        me.nedg[n] = toarr((x[12], x[15]))
        me.ncsp[n] = toarr((x[13], x[16]))
        me.fbsm[n] = toarr((x[14], x[18], x[17]))



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
        me.B[n] = (1. - n) * fac * z * fbav
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
        s[i] += "%14.0f %14.0f %14.0f " % (
            me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
        s[i] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
            me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
            me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
            me.B[n])
        if me.haserr:
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
      s[n] += "%14.0f %14.0f %14.0f " % (
          me.nedg[n][0], me.ncsp[n][0], me.fbsm[n][0])
      s[n] += "%18.14f %16.14f %16.14f %+17.14f %+20.14e " % (
          me.nedg[n][1] / me.nedg[n][0], me.ncsp[n][1] / me.ncsp[n][0],
          me.fbsm[n][2] / me.fbsm[n][0], me.fbsm[n][1] / me.fbsm[n][0],
          me.B[n])
      if me.haserr:
        s[n] += "%9.2e %9.2e %9.2e " % (
            me.eB[n], me.eZZ[n], me.eZn[n])
      s[n] = s[n].rstrip() + "\n"
      me.tot += me.hist[i]
    s[0] = "# %d %d V3 %d 1 %s\n" % (
        me.D, me.nmax, me.nmin, me.nedxmax)
    open(fn, "w").writelines(s)
    print "saved Zr file %s, tot %g" % (fn, me.tot)



def sumdat(fninp, fnls, fnout = None, fnZr = None):
  ''' sum over data in fninp + fnls, output to fnZrr '''
  gc = GC(fninp)
  m = len(fnls) + 1
  gcls = [None] * m
  fnls = [fninp,] + fnls
  # sum over all file lists
  for i in range(m):
    fn = fnls[i]
    gcls[i] = GC(fn)
    gcls[i].recompute()
    if i > 0: gc.absorb(gcls[i])
  gc.recompute()

  # compute errors
  Bls  = [None] * m
  Znls = [None] * m
  ZZls = [None] * m
  for i in range(m):
    Bls[i]  = gcls[i].B
    Znls[i] = gcls[i].Zn
    ZZls[i] = gcls[i].ZZ
  gc.eB  = errarr( Bls )
  gc.eZn = errarr( Znls )
  gc.eZZ = errarr( ZZls )
  gc.haserr = 1

  if not fnout: fnout = fninp + "a"
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



def getfninp():
  ''' get the input file '''
  if len(sys.argv) > 1:
    return sys.argv[1]
  ls = glob.glob("ZrrD*n*.dat")
  if len(ls): return ls[0]
  ls = glob.glob("ZrhD*n*.dat")
  if len(ls): return ls[0]



def getslaves(fninp):
  ''' get a list file names of non-master nodes '''
  ls = []
  i = 1
  while 1:
    fn = "%s%d" % (fninp, i)
    if not os.path.exists(fn):
      break
    ls += [fn,]
    i += 1
  return ls



# main function starts here
if __name__ == "__main__":
  fninp = getfninp()
  fnls = getslaves(fninp)
  sumdat(fninp, fnls)

