#!/usr/bin/env python

import os, sys, glob, re, getopt
from math import *



maxdim = 30



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -D:     maximal dimension
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "D:",
        [ "dim=",
          "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, maxdim

  for o, a in opts:
    if o in ("-D", "--dim"):
      maxdim = int(a)
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()



def getprec(s):
  if s == "mpq": return 100000
  elif s == "f128": return 113
  elif s == "ldbl": return 64
  else: return 53



def doit(dim):
  # 1. adjust the input list
  fns = glob.glob("gaussfn*mpf.dat")
 
  # 2. collect the available orders
  ls = []
  ns = []
  for fn in fns:
    # letter E means a partial input
    if "E" in fn: continue
    m = re.match("gaussfn([0-9]+)(.*)mpf.dat", fn)
    if not m:
      print fn
      raise Exception
    n = int( m.group(1) )
    ns += [n,]
    ls += [(n, getprec(m.group(2)), fn),]
  ns = sorted(list(set(ns)))
  
  # 3. handle orders
  output = "# %s\n" % (ns[-1])
  for n in ns:
    ls1 = [x for x in ls if x[0] == n]
    # choose the file with the highest precision
    n1, prec1, fn1 = sorted(ls1, key = lambda x: x[1])[-1]
    for s in open(fn1).readlines():
      if s.startswith("#"): continue
      d2, svir = s.strip().split()
      if int(d2) == dim:
        output += "%4s %s\n" % (n, svir)
        #print dim, n, float(svir)
        #raw_input()
        break
  
  # 4. save the output file
  fn = "gaussfD%dmpf.dat" % dim
  print "D %4s: writing %s" % (dim, fn)
  open(fn, "w").write(output)
 


if __name__ == "__main__":
  doargs()
  for dim in range(1, maxdim+1):
    doit(dim)

