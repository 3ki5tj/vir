#!/usr/bin/env python

''' write readable format for virXx.txt by
    convert it to BnPY.dat '''

import os, sys, re, glob



def conv(fninp, fnout):
  src = open(fninp).read().strip()
  # remove line continuation
  src = src.replace("\\\n", "")
  k = src.rfind("{")
  if k <= 0:
    return -1
  src = src[k:].strip("{}")
  src = src.replace(" ", "").replace("\n", "")
  src = src.split(",")
  # remove the trailing `40. after each figure
  for i in range(len(src)):
    s = src[i]
    m = re.match(r"(.*)(`[0-9]+\.)(.*)", s)
    if m:
      src[i] = m.group(1)
      tail = m.group(3).strip()
      if len(tail) > 2:
        #  *^(-6) ==> e-6
        src[i] += "e" + tail[2:]
    src[i] = float(src[i])

  out = "#  n   Bn/B2^(n-1)\n"
  for i in range(len(src)):
    out += "%4s %+23.15e\n" % (i + 1, src[i])
  print "converting %s to %s" % (fninp, fnout)
  open(fnout, "w").write(out)



fnls = glob.glob("vir*.txt")
for fninp in fnls:
  m = re.match("vir(.*).txt", fninp)
  fnout = "BnPY" + m.group(1) + ".dat"
  conv(fninp, fnout)
