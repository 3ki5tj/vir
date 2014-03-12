#!/usr/bin/env python

import os, sys, re, getopt



n = 64
D = 50
d = 3
nsteps = 2000
width = 800
height = 800



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "n:d:D:1:w:h:",
        ["nsteps=", "order=", "indim=", "outdim=", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, n, D, d, nsteps, width, height

  for o, a in opts:
    if o in ("-n", "--order"):
      n = int(a)
    elif o in ("-D", "--indim"):
      D = int(a)
    elif o in ("-d", "--outdim"):
      d = int(a)
    elif o in ("-1", "--scan"):
      nsteps = int(a)
    elif o in ("-w",):
      width = int(a)
    elif o in ("-h",):
      height = int(a)
    elif o in ("--verbose=",):
      verbose = int(a)
    elif o in ("--help",):
      usage()


# main function starts here
if __name__ == "__main__":
  doargs()
  os.system("make draw mds")
  fnmds = "mdsD%sn%sd%s.pos" % (D, n, d)
  if not os.path.exists(fnmds):
    fnpos = "D%sn%s.pos" % (D, n)
    if not os.path.exists(fnpos):
      os.system("icc -DD=%s -DN=%s mcrat0.c && ./a.out -L0 --nstpos=1000000 -w100000000 --hash-mode=0" % (D, n))
    os.system("./mds %s -d %s -1 %s" % (fnpos, d, nsteps))
  os.system("./draw %s -w %s -h %s" % (fnmds, width, height))

