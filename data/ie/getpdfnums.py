#!/usr/bin/env python




''' get numbers from text pasted from PDF file '''



from math import *



def geterr(x):
    ''' return the error from the significant digits '''
    x = fabs(x)
    if x == 0:
        return 0
    n = ceil(log10(x) + 1)
    while 1:
        eps = 10**n
        rem = fmod(x/eps + 0.5, 1.0) - 0.5
        if fabs(rem) < 0.00001:
            return eps
        n -= 1



def getpdfnums(fn):
    s = open(fn).readlines()
    arr = []
    err = []
    for i in range(len(s)):
        ln = s[i].rstrip()
        if ln == "":
            continue
        ln = ln.replace(" ", "")
        # replace the minus sign in utf-8
        ln = ln.replace("\xe2\x88\x92", "-")

        # skip the number in the parentheses (the error)
        p = ln.find("(")
        if p >= 0:
            q = ln.find(")")
            if q < 0:
                print ln
                raise Exception
            ln = ln[:p] + ln[q+1:]

        # find the times sign
        p = ln.find("\xc3\x97")
        if p >= 0:
            x = float(ln[:p])
            xp = ln[p+2:]
            if not xp.startswith("10"):
                print "bad line", ln
                raise Exception
            xp = int(xp[2:])
            x *= 10.0 ** xp
        else:
            x = float(ln)
        arr += [ x ]
        err += [ geterr(x) ]
    return arr, err



def compare(arr1, arr2, err2):
    n = len(arr2)
    for i in range(n):
        x = arr1[i]
        y = arr2[i]
        err = err2[i]
        #print x, y, err
        if fabs(x - y) > err:
            print "Please verify:", x, y, err, fabs(x-y)
            raw_input()
    print "Proofing finished!"




