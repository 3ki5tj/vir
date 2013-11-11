#!/usr/bin/env python



''' assemble a file for the partition function '''



import re, os, sys, glob



def loadZ(fn):
  ''' load the partition function '''
  lines = open(fn).readlines()
  info = lines[0][1:].split()
  dim = int(info[0])
  nmax = int(info[1])
  version = info[2]
  if version.startswith("V"):
    version = int(info[2][1:])
  else: # old version 0 data
    return None
  s = "%s\t" % dim
  for i in range(1, nmax + 1):
    ss = lines[i].split()
    n = ss[0]
    Zn = ss[2]
    s += "%s\t" % Zn
  s = s.strip() + "\n"
  return s



def mkZ(fnout):
  ''' compile the partition functions at different dimensions '''
  src = ""
  for dim in range(2, 1000):
    Zf = glob.glob("ZrD%sr[1-9]n*.dat" % dim)
    if len(Zf) == 0:
      Zf = glob.glob("ZrD%sn*.dat" % dim)
      if len(Zf) == 0: break
    Zf = Zf[0]
    s = loadZ(Zf)
    if s == None: break
    print "absorbing %s" % Zf
    src += s
  if os.path.exists(fnout):
    src0 = open(fnout).read()
    if src0 == src:
      print "no need to update %s\n" % fnout
      return
  print "updating %s\n" % fnout
  open(fnout, "w").write(src)



if __name__ == "__main__":
  mkZ("Z.dat")
