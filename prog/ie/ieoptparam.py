#!/usr/bin/env python



''' optimize parameters '''



import os, sys, glob, re, getopt
from math import *



dim = 6
nmin, nmax = 4, 12 # target range
gaussf = 0
shift_min, shift_max, shift_del = None, None, None
shiftn_min, shiftn_max, shiftn_del = None, None, None
shiftinc_min, shiftinc_max, shiftinc_del = None, None, None
shiftl0_min, shiftl0_max, shiftl0_del = None, None, None

program = ""
options = "" # command-line options to be passed on to iegsl or ieodfftw


refdir = None # directory for the reference files
fnref = None # file for reference values
verbose = 0



def getminmaxdel(s, isint = False):
  if ":" not in s: # a single value
    if isint:
      return int(s), int(s), 1
    else:
      return float(s), float(s) + .5, 1.0

  arr = s.split(":")
  if len(arr) == 2:
    arr = [arr[0], 1, arr[1]] # the increment is 1
  elif len(arr) != 3:
    print "%s is not of the format min:del:max" % s

  if isint:
    return int(arr[0]), int(arr[2]), int(arr[1])
  else:
    try:
      return float(arr[0]), float(arr[2]), float(arr[1])
    except:
      print "cannot handle %s, %s" % (s, arr)
      raise Exception



def usage():
  """ print usage and die """
  print sys.argv[0], "[Options]"
  print """
  Extrapolate virial coefficients

  OPTIONS:
   -D:        dimension
   -n:        minimal order
   -N:        maximal order
   -G:        Gaussian model
   -c:        shift, in min:del:max format
   -C:        shift, in terms of n, instead of n - 2, in min:del:max format
   -d:        shift increment, in min:del:max format
   -L:        shift minimal order l0, in min:max format
   -P:        program to run, e.g., iegsl_l
   -O:        command-line options to be passed to iegsl or ieodfftw
   --ref:     reference file
   --dir:     reference directory
  """
  exit(1)



def doargs():
  ''' Handle common parameters from command line options '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "D:n:N:GC:c:d:L:P:O:v",
        [ "dim=", "nmin=", "nmax=",
          "ref=", "dir=",
          "shift=", "shiftn=", "shiftinc=", "shiftl0=",
          "prog=", "program=", "options=", "option=", "opt=",
          "verbose=", "help", ])
  except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()

  global dim, nmin, nmax, gaussf, verbose
  global program, options
  global fnref, refdir
  global shift_min, shift_max, shift_del
  global shiftn_min, shiftn_max, shiftn_del
  global shiftinc_min, shiftinc_max, shiftinc_del
  global shiftl0_min, shiftl0_max, shiftl0_del

  for o, a in opts:
    if o in ("-D", "--dim"):
      dim = int(a)
    elif o in ("-n", "--nmin"):
      nmin = int(a)
    elif o in ("-N", "--nmax"):
      nmax = int(a)
    elif o in ("-G",):
      gaussf = 1
    elif o in ("--ref",):
      fnref = a
    elif o in ("--dir",):
      refdir = a
    elif o in ("-c", "--shift"):
      shift_min, shift_max, shift_del = getminmaxdel(a)
    elif o in ("-C", "--shiftn"):
      shiftn_min, shiftn_max, shiftn_del = getminmaxdel(a)
    elif o in ("-d", "--shiftinc"):
      shiftinc_min, shiftinc_max, shiftinc_del = getminmaxdel(a)
    elif o in ("-L", "--shiftl0"):
      shiftl0_min, shiftl0_max, shiftl0_del = getminmaxdel(a, True)
    elif o in ("-P", "--prog", "--program"):
      program = a
    elif o in ("-O", "--options", "--option", "--opt"):
      options = a
      if options.startswith('"'): options = options[1:]
      if options.endswith('"'): options = options[:-1]
      options = options.strip()
    elif o in ("-v",):
      verbose += 1
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("--help",):
      usage()

  if shiftl0_min == None:
    if gaussf:
      shiftl0_min, shiftl0_max, shiftl0_del = 3, 6, 1
    else:
      shiftl0_min, shiftl0_max, shiftl0_del = 2, 6, 1

  if shift_min == None:
    shift_min, shift_max, shift_del = 0.0, 0.5, 1.0

  if shiftinc_min == None:
    shiftinc_min, shiftinc_max, shiftinc_del = 0.0, 0.5, 1.0

  if not program:
    program = "./iegsl" if dim % 2 == 0 else "./ieodfftw_q -M8192"
  else:
    if not program.startswith("."): program = os.path.join(".", program)

  options = (options.rstrip() + " --corr").strip()

  print "dim %d, n [%d, %d], Gaussian %s, program %s, options %s" % (
      dim, nmin, nmax, gaussf, program, options)
  print "shift %s:%s:%s" % (shift_min, shift_del, shift_max)
  print "shiftinc %s:%s:%s" % (shiftinc_min, shiftinc_del, shiftinc_max)
  print "shiftl0 %s:%s:%s" % (shiftl0_min, shiftl0_del, shiftl0_max)
  #print "Enter to start",
  #raw_input()



def niceprint(x, err, dim = 0, n = 0, cnt = "", lsprec = ""):
  try:
    import scifmt
    print "%4d|%4d: %30s %+20.10e %9.2e %s %s" % (
        dim, n, scifmt.scifmt(x, err).text(errmax = errmax),
        x, err, cnt, lsprec)
  except Exception:
    print "%4d %24.14e %e %s %s" % (n, Bn, err, cnt, lsprec)



def loadref(dim, nmax, gaussf, fnref):
  if fnref == None:
    # guess the default
    if gaussf: # gaussian model
      fnref = "gaussfD%smpf.dat" % dim
      if refdir:
        fnref = os.path.join(refdir, fnref)
      if not os.path.isfile(fnref):
        print "no file %s for Gaussian model" % fnref
        raise Exception
    else:
      pat = "BnD%sn*.dat" % dim
      if refdir:
        pat = os.path.join(refdir, pat)
      fns = glob.glob(pat)
      if not fns:
        print "no reference for hard-sphere dim %d" % (dim)
        raise Exception
      nn = 1000
      for fn in fns:
        m = re.search(r"BnD%sn([0-9]+)\.dat" % dim, fn)
        if not m: continue
        n1 = int( m.group(1) )
        if n1 < nmax: continue
        if not fnref or n1 < nn: # choose a file with smaller nmax
          nn, fnref = n1, fn
      if not fnref:
        print "no suitable reference file for the hard-sphere model, only %s" % fns
        raise Exception

  refs = [0] * (nmax + 1)
  # load the content of the reference file
  for s in open(fnref).readlines():
    if s.strip() == "" or s.startswith("#"): continue
    a = s.strip().split()
    n1 = int(a[0])
    if n1 >= 1 and n1 <= nmax:
      refs[n1] = float(a[1])

  return refs



def testshift(dim, nmin, nmax, gaussf,
    shifttype, shift, shiftinc, shiftl0, refs, verbose = False):
  # return the maximal relative error in the range
  #print shift, shiftinc, shiftl0

  # remove old files
  pat = "BnPYcD%sn%s*.dat" % (dim, nmax)
  if gaussf: pat = "GF" + pat
  if dim % 2 == 0: pat = "h" + pat
  os.system("rm -f " + pat)

  # run the command
  sys = "-G" if gaussf else ""
  cmd = "%s %s %s -D%s -n%s -%s%s -d%s -L%s" % (
      program, options, sys, dim, nmax,
      shifttype, shift, shiftinc, shiftl0)
  if not verbose:
    cmd += " 2>ieoptparamD%sn%s.err >ieoptparamD%sn%s.out" % (
        dim, nmax, dim, nmax)
  if verbose: print cmd
  os.system(cmd)

  # load the output
  try:
    fnout = glob.glob(pat)[0]
  except:
    print "cannot find %s" % pat
    raise Exception
  errmax = errsum = cnt = 0
  who = None
  for ln in open(fnout).readlines():
    if ln.strip() == "" or ln.startswith("#"): continue
    a = ln.strip().split()
    n1 = int(a[0])
    if n1 < nmin or n1 > nmax: continue
    x = float(a[3])
    ref = refs[n1]
    # compute the relative error
    err = fabs(x / ref - 1)
    if x * ref <= 0:
      # penalize wrong signs
      err = fabs(x / ref - 1) * 100
    else:
      err = max(fabs(log(x / ref)), fabs(x/ref - 1))
    ## we let higher-order coefficients carry larger weight
    ## in computing the average error
    #wt = n1 * n1
    wt = 1
    errsum += err * wt
    cnt += wt
    if err > errmax:
      errmax = err
      who = n1
    if verbose: print x, ref, err, errmax
  if not who:
    print "bad result:", cmd, "\n", open(fnout).read()
    raise Exception
  #raw_input()
  errav = errsum/cnt
  # based on the maximal error
  return errmax, errav, who
  # based on the average error
  #return errav, errmax, who



def optshift(dim, nmin, nmax, gaussf, fnref):
  # optimize the shift parameters

  refs = loadref(dim, nmax, gaussf, fnref)
  print "reference values", refs
  errmin, who = 1e10, None
  erravmin = 1e10

  shift_best, shiftinc_best, shiftl0_best = None, None, None

  opt0 = program + " --corr" + (" -G" if gaussf else "") + (" -D%s" % dim)

  if shiftn_min != None:
    shifttype = "C"
  else:
    shifttype = "c"
  strsh = ""

  for shiftl0 in range(shiftl0_min, shiftl0_max + 1):

    if shifttype == "C":
      shift = shiftn_min
    else:
      shift = shift_min

    while 1:

      shiftinc = shiftinc_min
      while shiftinc < shiftinc_max:
        err, errav, who = testshift(dim, nmin, nmax, gaussf,
            shifttype, shift, shiftinc, shiftl0, refs)
        # if the minimal error are the same, compare the average error
        if err < errmin or fabs(err/errmin - 1) < 0.001 and errav < erravmin:
          shift_best, shiftinc_best, shiftl0_best = shift, shiftinc, shiftl0
          if shifttype == "C":
            strsh = " (c %s)," % (shift + shiftinc*2)
          else:
            strsh = " (c' %s)," % (shift - shiftinc*2)
          print "* best parameters: %s -%s%s -d%s -L%s,%s error %s%% (n %s, %s%%)" % (
              opt0, shifttype, shift, shiftinc, shiftl0, strsh, err*100, who, errav*100)
          errmin, erravmin, whomin = err, errav, who
        elif verbose:
          print "parameters: %s -%s%s -d%s -L%s error %s%% (n %s, %s%%)" % (
            opt0, shifttype, shift, shiftinc, shiftl0, err*100, who, errav*100)
        shiftinc += shiftinc_del

      if shifttype == "C":
        shift += shiftn_del
        if shift > shiftn_max: break
      else:
        shift += shift_del
        if shift > shift_max: break

  testshift(dim, nmin, nmax, gaussf,
      shifttype, shift_best, shiftinc_best, shiftl0_best, refs, True)
  print "* best parameters: %s -%s%s -d%s -L%s,%s error %s%% (n %s, %s%%)" % (
      opt0,
      shifttype, shift_best, shiftinc_best, shiftl0_best, strsh,
      errmin*100, whomin, erravmin*100)


if __name__ == "__main__":
  doargs()
  optshift(dim, nmin, nmax, gaussf, fnref)


