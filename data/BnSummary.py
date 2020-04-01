#!/usr/bin/env python

import glob, scifmt, re
from math import *


def add_Bn(fn, tab):
    for s in open(fn).readlines():
        if s[0] == '#':
            continue
        x = s.strip().split('\t')
        n = int(x[0])
        Bn = float(x[1])
        err = float(x[2])
        dat = (Bn, err)
        if not n in tab:
            tab[n] = [dat,]
        else:
            relerr = err/fabs(Bn)
            if dat in tab[n] and relerr > 1e-15:
                print ("{}:{} already exists".format(fn, dat))
            else:
                tab[n] += [dat,]


def summarize_Bn(dim, tab):
    nmax = max( [n for n in tab] )
    out = [1]*(1 + nmax)
    for n in tab:
        x = tab[n]
        m = len(x) # number of virial coeffient
        y = x[0][0]
        sig = x[0][1]
        if x[0][1] > 1e-300 and m > 1:
            ytot = 0
            wtot = 0
            for j in range(m):
                y = x[j][0]
                sig = x[j][1]
                w = 1.0/(sig*sig)
                ytot += y * w
                wtot += w
                if n == 7 and dim == 100:
                    print "D100n7", j, sig, w, sqrt((j+1)/wtot)
            y = ytot / wtot
            sig = (1.0 / wtot) ** 0.5
        out[n] = (y, sig)
        #print n, y, sig, m
        #raw_input()
    
    s = "{}\t1\t1".format(dim)
    for n in range(3, len(out)):
        y, sig = out[n]
        s += "\t" + scifmt.scifmt(y, sig).text()
    s += '\n'

    return s


def handle_dim(d):
    fns = glob.glob("BnD{}n*.dat".format(d))
    if not fns: return ""

    # collect virial coefficients from different files
    tab = {}
    for fn in fns:
        add_Bn(fn, tab)

    # summarize virial coefficients
    s = summarize_Bn(d, tab)

    return s


# output 
s = "D"
for n in range(1, 64):
    s += "\tn=" + str(n)
s += '\n'

for d in range(1, 101):
    s += handle_dim(d)

fnout = "BnSummary.txt"
open(fnout, "w").write(s)
print ("saved output to " + fnout)

fncsv = fnout.replace(".txt", ".csv")
s_csv = s.replace('\t', ',')
open(fncsv, "w").write(s_csv)
print ("saved output to " + fncsv)

