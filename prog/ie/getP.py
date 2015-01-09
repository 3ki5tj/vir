#!/usr/bin/env python

import sys, os, re
from math import *



def getP(vir, den):
  return sum(vir[i]*den**i for i in range(len(vir)))



def getmu(vir, den):
  return sum((i+1.)/i*vir[i+1]*den**i for i in range(1, len(vir)-1))



def getvir(fn, col = 1):
  n = len(open(fn).readlines())
  vir = [0]*(n + 4)
  for ln in open(fn).readlines():
    if ln.startswith("#"): continue
    arr = ln.strip().split()
    order = int(arr[0])
    x = float(arr[col])
    vir[order] = x
  vir[1] = 1
  vir[2] = 1
  return vir



def normalize(vir, B2 = None):
  if B2 == None: B2 = 2*pi/3
  #print "B2 %g" % B2
  for i in range(2, len(vir)):
    vir[i] *= pow(B2, i - 1)



def demo():
  # Hutchinson-Conkie, 3D hard-sphere
  vir = [0, 1, 2.0943951024e+00,
   2.7415568e+00,
   2.5939709e+00,
   2.1529229e+00,
   1.6182668e+00,
   1.1381708e+00,
   8.2357838e-01,
   4.9375454e-01,
   3.5481524e-01,
   1.9112390e-01,
   6.3862960e-02,
   1.7663569e-01,
  -9.7172953e-02,
   1.9196284e-01,
  -1.0236569e-02,
  -2.0923830e-01,
   6.0471239e-01,
  -9.8117704e-01,
   1.0758527e+00,
  -4.3738719e-01,
  -2.0923830e-01,
  # 6.0471239e-01,
  #-9.8117704e-01,
  # 1.0758527e+00,
  #-4.3738719e-01,
  #-1.3100728e+00,
  # 4.2300029e+00,
  #-7.5002981e+00,
  ]

  for i in range(1, 15):
    den = 0.05*i
    rho = den * sqrt(2)
    print "%.2f %s" % (den, getP(vir, rho)/rho - 1)



def doit(fn, col, n1 = -1, n2 = -1):
  vir = getvir(fn, col)
  normalize(vir)
  text = "#  rho  beta*P/rho-1     beta*P            beta*mu\n"
  for i in range(1, 1002):
    rho = 0.001*i
    if n1 > 0 and n2 > 0:
      # here, we average two truncated series
      # at orders n1 and n2
      vir1 = vir[:n1+1]
      vir2 = vir[:n2+1]
      pres = (getP(vir1, rho) + getP(vir2, rho))/2
      mu = (getmu(vir1, rho) + getmu(vir2, rho))/2
    else:
      pres = getP(vir, rho)
      mu = getmu(vir, rho)
    line = "%.4f %-16s %-16s %-16s" % (rho, pres/rho - 1, pres, mu)
    text += line.rstrip() + "\n"
  print text



if __name__ == "__main__":
  if len(sys.argv) > 1:
    fnin = sys.argv[1]
    # if the file contains multiple columns
    # for different types of the virial coefficients
    # we can use the second command-line argument
    # to determine the column id
    if len(sys.argv) > 2:
      col = int(sys.argv[2])
    else:
      col = 1
    # here we want to do two truncated series
    # at n1 and n2
    if len(sys.argv) >= 5:
      n1 = int(sys.argv[3])
      n2 = int(sys.argv[4])
    else:
      n1 = n2 = -1
    doit(fnin, col, n1, n2)
  else:
    demo()

