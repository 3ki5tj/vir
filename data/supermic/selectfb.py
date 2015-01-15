#!/usr/bin/env python


import os, sys


def select(fn1, fn2):
  ''' select the larger file of fn1 and fn2 '''
  if not os.path.exists(fn1) or not os.path.exists(fn2):
    return
  sz1 = os.path.getsize(fn1)
  sz2 = os.path.getsize(fn2)
  if sz1 < sz2:
    print "replacing %s by %s" % (fn1, fn2)
    os.rename(fn2, fn1)



def selectfb():
  select("fb.bdb", "fb.bdb.bak")
  select("mic0/fb.bdb", "mic0/fb.bdb.bak")
  select("mic1/fb.bdb", "mic1/fb.bdb.bak")



if __name__ == "__main__":
  selectfb()
