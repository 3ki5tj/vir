#!/usr/bin/env python



import glob, os, sys
import scifmt
import getpdfnums



def loadvir(fn, B = None, err = None):
    if not B:   B   = [0]*132
    if not err: err = [0]*132
    for s in open(fn).readlines():
        if s.startswith("#"): continue
        a = s.split()
        n = int(a[0])
        Bn = float(a[1])
        errn = float(a[2])
        if B[n] == 0 or err[n] > errn:
            B[n] = Bn
            err[n] = errn
    return B, err



def mkrow(dim):
    # integral equation
    fn = glob.glob("../xBnPYcD%sn*.dat" % dim)[0]
    B, err = loadvir(fn)

    # Mayer sampling
    fn2ls = glob.glob("../../BnD%sn*.dat" % dim)
    B2, err2 = loadvir(fn2ls[0])
    for fn2 in fn2ls[1:]:
        B2, err2 = loadvir(fn2, B2, err2)

    arr = []
    for n in range(32, 129, 32):
        if B2[n] != 0:
            arr += [ B[n], B2[n] ]
        else:
            arr += [ B[n] ]

    return arr



def checkTableIII():
    arrpdf, errpdf = getpdfnums.getpdfnums("TableIIIpdf.txt")
    arr = []
    for n in range(10, 31):
        arr += mkrow(n)

    getpdfnums.compare(arr, arrpdf, errpdf)



if __name__ == "__main__":
    checkTableIII()
