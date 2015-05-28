#!/usr/bin/env python



import glob, os, sys, re
import scifmt
import gaussfextra
import getpdfnums



def loadvir(fn):
    B   = [0]*101
    err = [0]*101

    if fn.startswith("GFBn") or fn.startswith("hGFBn"):
        # DSC
        colB = 3
    else:
        # exact
        colB = 1

    for s in open(fn).readlines():
        if s.startswith("#"): continue
        a = s.split()
        n = int(a[0])
        B[n] = float(a[ colB ])
        err[n] = abs(B[n]) * 1e-20
    return B, err



def getref(dim, nmin, nmax):
    ''' get reference data '''
  
    # DSC
    a = glob.glob("*GFBnPYcD%dn*.dat" % dim)
    if len(a[0]) < len(a[1]):
        fn1, fn2 = a[0], a[1]
    else:
        fn1, fn2 = a[1], a[0]
    B1, err1 = loadvir(fn1)
    B2, err2 = loadvir(fn2)
  
    arr = []
    for n in range(nmin, nmax + 1):
        # exact result
        Bx, errx = gaussfextra.extrapolate(dim, n)
        arr += [B1[n], B2[n], Bx]
    return arr



def checkTableV():
    arrpdf, errpdf = getpdfnums.getpdfnums("TableVpdf.txt")

    arr = []
    arr += getref(8,  13, 16)
    arr += getref(9,  13, 16)
    arr += getref(9,  17, 20)
    arr += getref(10, 13, 16)
    arr += getref(10, 17, 20)

    getpdfnums.compare(arr, arrpdf, errpdf)



if __name__ == "__main__":
    checkTableV()
