#!/usr/bin/env python



'''
automate the process of visualizing high-dimensional clusters
calls mcrat0.c for sampling,
mds.c for multidimensional scaling
and draw.c for the visualization

The basic usage is
  python mdsdraw.py -D 50 -n 64 -d 3
This command computes 50 dimensional virial coefficients
of the 64th order.  The resulting clusters are visualized
in three dimensions.  The command took several minutes to finish.

Cf. README.mds for details.
'''



import os, sys, re, getopt



n = 64
D = 50
d = 3
nsteps = 2000
width = 800
height = 800
forceMDS = 0
ampMDS = 1.0
dtMDS = 0.01


def usage():
  print "%s [Options] your.pos"



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "n:d:D:1:w:h:a:t:",
        ["nsteps=", "order=", "indim=", "outdim=",
         "MDS", "verbose=", "help",])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global verbose, n, D, d, nsteps, width, height
  global ampMDS, dtMDS, forceMDS

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
    elif o in ("-a",):
      ampMDS = float(a)
    elif o in ("-t",):
      dtMDS = float(a)
    elif o in ("--MDS",):
      forceMDS = 1
    elif o in ("--verbose=",):
      verbose = int(a)
    elif o in ("--help",):
      usage()


# main function starts here
if __name__ == "__main__":
  doargs()
  os.system("make draw mds")
  fnmds = "mdsD%sn%sd%s.pos" % (D, n, d)
  if not os.path.exists(fnmds) or forceMDS:
    fnpos = "D%sn%s.pos" % (D, n)
    if not os.path.exists(fnpos):
      cmd = ("icc -DD=%s -DN=%s mcrat0.c " % (D, n)
          + "&& ./a.out -L0 --nstpos=1000000 -w100000000 --hash-mode=0")
      os.system( cmd )
    os.system("./mds %s -d %s -1 %s -a %s -t %s" % (fnpos, d, nsteps, ampMDS, dtMDS))
  os.system("./draw %s -w %s -h %s" % (fnmds, width, height))

