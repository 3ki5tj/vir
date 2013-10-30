#!/usr/bin/env python



def de_bruijn(k, n):
  """
  De Bruijn sequence for alphabet size k
  and subsequences of length n.
  this code is copied from wikipedia
  """
  a = [0] * k * n
  sequence = []
  def db(t, p):
    if t > n:
      if n % p == 0:
        for j in range(1, p + 1):
          sequence.append(a[j])
    else:
      a[t] = a[t - p]
      db(t + 1, p)
      for j in range(a[t - p] + 1, k):
        a[t] = j
        db(t + 1, t)
  db(1, 1)
  return sequence


def bits2int(s):
  i = 0
  for b in s:
    i = i * 2 + b
  return i



def getarr(s, n):
  arr = [0] * (2**n)
  for i in range(2**n):
    ss = s[i:i+n]
    if len(ss) < n:
      ss = ss + [0] * (n - len(ss))
    q = bits2int(ss)
    arr[q] = i
    print "%4d %3d %s" % (i, q, ss)
  return arr


# verify the magic number works
def verify(n, magic, arr):
  pow2n = 2**n # n = 6, pow2n 64
  pow2to2n = 2**pow2n  # 2^64
  mask = pow2to2n - 1
  print "0x%x" % mask
  for i in range(pow2n): # 0..63
    b = 2**i
    x = magic * b
    x = x & mask
    low = x >> (pow2n - n) # 58
    print "%2d %3d %2d" % (i, low, arr[low])
    if arr[low] != i:
      raise Exception

n = 6
s = de_bruijn(2, n)
magic = bits2int(s)
arr = getarr(s, n)
verify(n, magic, arr)
print "0x%x %s" % (magic, arr)

