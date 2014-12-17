#!/usr/bin/env python



'''
clean up snapshot files for y(r), snapshot_yr*.dat,
produced by ievir.c (or iegsl.c) with -DDHT

Since the discrete Hankel transform (DHT) for even dimensions
  is slow, snapshots are used for restartable runs.
The snapshot files for different correlation functions,
  e.g., c(r) and y(r), however, may not be synchronized.
By default, for each order n, a file for y(r) is produced,
  although usually only those for the last two orders are needed
  and this script clears these files for low orders.
'''



import sys, os, glob, re, getopt



nkeep = 2 # number of recent items to keep
pattern = "snapshot*_yrl*.dat"



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options] [pattern]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -k:     number of recent items to keep
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "k:v",
        [ "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global nkeep, pattern

  for o, a in opts:
    if o in ("-k",):
      nkeep = int(a)
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help",):
      usage()

  if len(args) > 0:
    pattern = args[0]
  #print pattern


def foo():
  fns = glob.glob(pattern)

  # 1. collect stems
  stems = list(set([fn.split("_")[0] for fn in fns]))
  #print stems

  # 2. clear yr
  for stem in stems:
    fnyrs = [fn for fn in fns if fn.split("_")[0] == stem]
    #print stem, fnyrs
    ls = []
    for fnyr in fnyrs:
      m = re.search(r"yrl([0-9]+)\.dat", fnyr)
      if not m: raise Exception
      ls += [(int(m.group(1)), fnyr),]
    ls = sorted(ls, key = lambda x: x[0])
    #print ls
    if len(ls) >= nkeep:
      for id, fn in ls[:-nkeep]:
        print "removing", fn
        os.remove(fn)


if __name__ == "__main__":
  doargs()
  foo()
