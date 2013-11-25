#!/usr/bin/env python



''' assemble a file for the partition function '''



import re, os, sys, glob



def loadZ(fn, Zarr = []):
  ''' load the partition function '''
  lines = open(fn).readlines()
  info = lines[0].split()[1:]
  dim = int(info[0])
  nmax = int(info[1])
  version = info[2]
  if version.startswith("V"):
    version = int(info[2][1:])
  else: # old version 0 data
    return None
  lcnt = min(nmax + 1, len(lines))
  for i in range(1, lcnt):
    ss = lines[i].split()
    n = int(ss[0])
    Zn = ss[2]
    if n >= len(Zarr): # extend the array
      Zarr += [0] * (n + 1 - len(Zarr))
    Zarr[n] = Zn
  return Zarr



def getnmax(fn):
  ''' deduce nmax '''
  m = re.search(r"n([0-9]+)", fn)
  if not m:
    print "cannot find nmax in %s" % fn
    raise Exception
  return int( m.group(1) )



def mysort(arr):
  ''' sort file names by the nstfb '''
  n = len(arr)
  for i in range(n): # settle position of i
    imax = i
    nmaxi = getnmax(arr[imax])
    for j in range(i+1, n):
      nmaxj = getnmax(arr[j])
      if nmaxj > nmaxi:
        imax = j
        nmaxi = nmaxj
    if imax != i: # swap i and imax
      arr[imax], arr[i] = arr[i], arr[imax]
  return arr



def Zarr2str(dim, arr):
  ''' array to string '''
  if len(arr) <= 1:
    return ""
  s = "%s\t" % dim
  for i in range(1, len(arr)):
    s += "%s\t" % arr[i]
  return s.strip() + '\n'



def dodim(dim):
  Zfls = glob.glob("ZrD%sr[1-9]n*.data" % dim)
  Zfls += glob.glob("ZrhD%sn*.data" % dim)
  Zfls = mysort(Zfls)
  src = ""
  Zarr = [0,]
  for Zf in Zfls:
    Zarr = loadZ(Zf, Zarr)
    print "absorbing %s" % Zf
  return Zarr2str(dim, Zarr)



def mkZ(fnout):
  ''' compile the partition functions at different dimensions '''
  src = ""
  for dim in range(2, 1000):
    src += dodim(dim)

  # write the output fnout
  if os.path.exists(fnout):
    src0 = open(fnout).read()
    if src0 == src:
      print "no need to update %s\n" % fnout
      return
  print "updating %s\n" % fnout
  open(fnout, "w").write(src)



if __name__ == "__main__":
  mkZ("Z.dat")
