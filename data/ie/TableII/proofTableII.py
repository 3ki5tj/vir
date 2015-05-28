#!/usr/bin/env python



import glob, os, sys
import scifmt
import getpdfnums



def loadvir(fn):
  Bc = [0]*15
  Bv = [0]*15
  By = [0]*15
  for s in open(fn).readlines():
    if s.startswith("#"): continue
    a = s.split()
    n = int(a[0])
    Bc[n] = float(a[1])
    Bv[n] = float(a[2])
    if fn.startswith("x"):
      By[n] = float(a[len(a)/2])
    else:
      By[n - 1] = float(a[-1])
  return Bc, Bv, By



def mkrow(name, tag, doy):
    fn = glob.glob("xBn" + tag + "D3*.dat")
    if len(fn) > 0:
        fn = fn[0]
    else:
        fn = glob.glob("Bn" + tag + "D3*.dat")[0]
    Bc, Bv, By = loadvir(fn)

    arr = []
    for n in range(4, 13):
        arr += [ Bc[n], Bv[n] ]
        if doy:
            arr += [ By[n] ]
    return arr



def checkTableII():
    arrpdf, errpdf = getpdfnums.getpdfnums("TableIIpdf.txt")

    arr = []
    arr += mkrow("YBG",          "YBG",  True)
    arr += mkrow("Kirkwood",     "K",    True)
    arr += mkrow("PY",           "PY",   True)
    arr += mkrow("HNC",          "HNC",  True)
    arr += mkrow("Hurst",        "H",    False)
    arr += mkrow("Rowlinson 1",  "R",    False)
    arr += mkrow("Rowlinson 2",  "IR",   False)
    arr += mkrow("HC",           "HC",   False)
    arr += mkrow("MS",           "MS",   False)
    arr += mkrow("BPGG",         "BPGG", False)
    arr += mkrow("Verlet",       "V",    False)
    arr += mkrow("MP",           "MP",   False)
    arr += mkrow("RY",           "RY",   False)
    arr += mkrow("Quadratic",    "SQR",  False)

    getpdfnums.compare(arr, arrpdf, errpdf)


if __name__ == "__main__":
    checkTableII()
