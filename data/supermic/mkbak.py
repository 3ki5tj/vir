#!/usr/bin/env python



import os, sys, glob
import selectfb



def mklog():
  logs = glob.glob("log/log*")
  if len(logs):
    ilog = 1 + max(int(d[7:]) for d in logs)
  else:
    ilog = 1
  fn = glob.glob("D*.out")
  if len(fn) > 0:
    fn = fn[0]
    fnlog = "log/log%s" % ilog
    print "moving %s to %s" % (fn, fnlog)
    os.system("mv %s %s" % (fn, fnlog))
  return ilog



def mkbak(ibak = -1):
  if ibak < 0 or os.path.exists("bak%s" % ibak):
    baks = glob.glob("bak*")
    if len(baks):
      ibak = 1 + max(int(d[3:]) for d in baks)
    else:
      ibak = 1
  dir = "bak%d" % ibak
  print "backing mr*.dat up to %s" % dir
  os.makedirs(dir + "/mic0")
  os.makedirs(dir + "/mic1")
  os.system("cp mr* " + dir)
  os.system("cp mic0/mr* " + dir + "/mic0")
  os.system("cp mic1/mr* " + dir + "/mic1")
  return ibak



if __name__ == "__main__":
  id = mklog()
  mkbak(id)
  selectfb.selectfb()
