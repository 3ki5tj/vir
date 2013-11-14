#!/usr/bin/env python

''' compute the equation of states from the virial coefficients
    run this script under data '''


drho = 0.001
nmax = 63

def getvirpy(fn):
  ''' obtain the virial coefficients from the Percus-Yevick results '''
  lines = open(fn).readlines()[1:]
  n = len(lines)
  vir = [1] * (n + 1)
  for i in range(n):
    order, bn = lines[i].strip().split()
    order = int(order)
    bn = float(bn)
    if order >= len(vir):
      print "%s contains large order %d virial coefficients" % (fn, order)
      raise Exception
    vir[order] = bn
  return vir



def getcompfac(vir, rho):
  ''' compute the compressibility factor beta*P/rho
      from the virial series for a particular rho '''
  nmax = len(vir)
  n = nmax - 1
  bPrho = 0
  while n >= 2:
    bPrho = rho * (vir[n] + bPrho)
    n -= 1
  return bPrho 



def geteos(vir, rhomax = 0.5, drho = 0.001):
  ''' get the equation of states from the virial coefficients
      note the rho is the density reduced by B2, i.e., rho* = rho*B2 '''
  nmax = int(rhomax/drho + .5) + 1
  bPrho = []
  for i in range(nmax):
    rho = i * drho
    bPrho += [ (rho, getcompfac(vir, rho)), ]
  return bPrho



def writeeos(eos, fn):
  ''' save equation of states '''
  nmax = len(eos)
  s = "#  rho*B2   beta*P/rho\n"
  for i in range(nmax):
    rho, vir = eos[i]
    s += "%8.6f %+22.14e\n" % (rho, vir)
  open(fn, "w").write(s)



def getpyeos(dim, cp, nmax = 100000, rhomax = 0.5, drho = 0.001):
  fninp = "pyhs/BnPY%d%s.dat" % (dim, cp)
  vir = getvirpy(fninp)
  if nmax >= len(vir):
    nmax = len(vir) - 1
  vir = vir[: nmax + 1]
  eos = geteos(vir, rhomax, drho)
  fnout = "eosPY%d%s.dat" % (dim, cp)
  writeeos(eos, fnout)
  return vir



def getvirZr(fnZr):
  ''' get the virial coefficients from a Zr file'''
  lines = open(fnZr).readlines()
  nmax = len(lines) - 1
  vir = [1] * (nmax + 1)
  for i in range(3, nmax + 1):
    s = lines[i].split()
    vir[i] = float(s[18])
  return vir

#print getvirZr("stampede/D13r1n64/ZrD13r1n64.data")



def getvirMR(dim):
  ''' get the virial coefficients from a directory like D2, D3, ...'''
  # get B3
  vir = [1] * 65
  return vir



def getsimeos(dim, fnZr, nmax = 10000, rhomax = 0.5, drho = 0.001):
  ''' get the equation of states from simulation data '''
  if fnZr:
    vir = getvirZr(fnZr)
  else:
    vir = getvirMR(dim)
  if nmax < len(vir) - 1:
    vir = vir[:nmax + 1]
  eos = geteos(vir, rhomax, drho)
  fnout = "eossim%d.dat" % dim
  writeeos(eos, fnout)
  return vir


def cmpeos3(dim, nmax = 10000, rhomax = 0.5, drho = 0.001):
  fnZr = "stampede/D%sr1n64/ZrD%sr1n64.data" % (dim, dim)
  virs = getsimeos(dim, fnZr, nmax = nmax, rhomax = rhomax)
  nmax = len(virs) - 1
  virc = getpyeos(dim, "c", nmax = nmax, rhomax = rhomax)
  virp = getpyeos(dim, "p", nmax = nmax, rhomax = rhomax)
  print len(virs), len(virc), len(virp)


cmpeos3(11)
