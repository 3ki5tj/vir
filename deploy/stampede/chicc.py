#!/usr/bin/env python

# change compiling options for icc

import os, shutil, re, sys




def mkbak(fn):
  i = 1
  while i < 100:
    fnbak = fn + ".bak" + str(i)
    if not os.path.exists(fnbak):
      shutil.copyfile(fn, fnbak)
      break
    i += 1


def chicc(fn, newopt):
  ''' change the icc command in the file `fn' '''
  if not os.path.exists(fn): return
  lines = open(fn).readlines()
  for i in range(len(lines)):
    s = lines[i]
    pos = s.find("icc")
    if pos < 0: continue

    while 1: # remove exisiting version
      pos1 = s.find(newopt)
      if pos1 < 0: break
      s = s[:pos1].rstrip() + " " + s[pos1 + len(newopt):].lstrip()

    pos = s.find(".icc")
    if pos < 0: continue
    s = s[:pos + 3] + " " + newopt + " " + s[pos + 3:].lstrip()
    lines[i] = s
    print fn, ": ", s,
    # back up
    mkbak(fn) ####
    open(fn, "w").writelines(lines)
    break


newopts = ["-ipo",]
for d in os.walk(os.curdir).next()[1]:
  for newopt in newopts:
    chicc( os.path.join(d, "foo.pbs"), newopt )
    chicc( os.path.join(d, "Makefile"), newopt )
