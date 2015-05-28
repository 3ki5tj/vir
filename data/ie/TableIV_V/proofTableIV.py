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



def getref(dim):
    ''' get reference data '''
  
    # exact result
    fn = "../../gaussf/gaussfD%smpf.dat" % dim
    Bx, errx = loadvir(fn)
  
    # DSC
    a = glob.glob("*GFBnPYcD%dn*.dat" % dim)
    if a:
        if len(a[0]) < len(a[1]):
            fn1, fn2 = a[0], a[1]
        else:
            fn1, fn2 = a[1], a[0]
        B1, err1 = loadvir(fn1)
        B2, err2 = loadvir(fn2)
        m = re.search(r"c([0-9\.]+)L(.*)ldbl", fn2)
        kappa = float( m.group(1) )
        lval = int( m.group(2) )
  
    arr = []
    for n in range(9, 13):
        if a:
            arr += [ B1[n], B2[n] ]
  
        # exact result
        arr += [ Bx[n] ]
  
    if a: # add a column for kappa
        arr += [ kappa, lval + 2 ]
  
    return arr



def checkTableIV():
    arrpdf, errpdf = getpdfnums.getpdfnums("TableIVpdf.txt")

    arr = []
    for d in range(1, 11):
        arr += getref(d)

    getpdfnums.compare(arr, arrpdf, errpdf)


if __name__ == "__main__":
  checkTableIV()
