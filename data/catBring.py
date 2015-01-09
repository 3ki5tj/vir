#!/usr/bin/env python

s1 = [s for s in open("Bring.dat").readlines() if s.strip() != ""]
s2 = [s for s in open("Bring2.dat").readlines() if s.strip() != ""]

i0 = 0  # offset from the main file
for i in range(len(s2)):
  dim1 = s1[i + i0].strip().split()[0]
  dim2 = s2[i].strip().split()[0]
  if dim1 != dim2:
    print "dimension mismatch %d vs %d" % (dim1, dim2)
    raise Exception
  s1[i + i0] = '\t'.join( s1[i + i0].strip().split()
                        + s2[i].strip().split()[1:] ) + '\n'
  #print s1[i + i0],

fnout = "Bring1.dat"
print "saved results to %s" % fnout
open(fnout, "w").writelines(s1)

