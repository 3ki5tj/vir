#!/usr/bin/env python



import glob, os, sys, re
import scifmt
import gaussfextra



def loadvir(fn):
  B   = [0]*101
  err = [0]*101

  if fn.startswith("GFBn") or fn.startswith("hGFBn"):
    # DSC
    colB = 3
  else:
    # exact
    colB = 1

  for s in open(fn).readlines():
    if s.startswith("#"): continue
    a = s.split()
    n = int(a[0])
    B[n] = float(a[ colB ])
    err[n] = abs(B[n]) * 1e-20
  return B, err



def getnum(x, err, errmax = 20):
  if x == 0:
    return "&nbsp;"

  it = 0
  while 1:
    s = scifmt.scifmt(x, err).html(errmax = errmax, xpmin = -2)
    i = s.find(".")
    j = s.find("(")
    jj = s.find(")")
    if s.startswith("&minus;"):
      k = 7
    elif s.startswith("+"):
      k = 1
    else:
      k = 0

    # comupte the number of significant digits
    ndigits = j - i - 1
    if s[k:i] != "0":
      ndigits += i - k
    if s[i+1] == "0":
      ndigits -= 1
      if i+2 < len(s) and s[i+2] == "0":
        ndigits -= 1

    # get the number of error digits
    nderr = jj - j - 1

    if ndigits < 10 + nderr:
      # remove the error estimation
      s = s[:j] + s[jj+1:]
      break
    else:
      err *= 10
      errmax = 10
    it += 1
  return s



def mkrow(dim, nmin, nmax):
  ''' write an HTML row '''

  srow = "<tr>\n"

  # DSC
  a = glob.glob("*GFBnPYcD%dn*.dat" % dim)
  print a, dim, nmin, nmax
  if len(a[0]) < len(a[1]):
    fn1, fn2 = a[0], a[1]
  else:
    fn1, fn2 = a[1], a[0]
  B1, err1 = loadvir(fn1)
  B2, err2 = loadvir(fn2)
  m = re.search(r"c([0-9\.]+)L(.*)ldbl", fn2)
  kappa = float( m.group(1) )
  lval = int( m.group(2) )

  for n in range(nmin, nmax + 1):
    srow += "<td>%d<sub>%d</sub></td>\n" % (dim, n)
    srow += "<td style='text-align: left;'>"

    # exact result
    Bx, errx = gaussfextra.extrapolate(dim, n)
    errx *= 2 # to be conservative
    errx = max(errx, abs(Bx) * 1e-9)

    # DSC, kappa = 0
    srow += getnum(B1[n], abs(B1[n] - Bx)) + "<br>"

    # DSC, kappa != 0
    srow += getnum(B2[n], abs(B2[n] - Bx)) + "<br>"

    srow += getnum(Bx, errx, errmax = 10)

    srow += "</td>\n"

  srow += "</tr>\n\n"
  return srow



def mkTableV():
  stab = "<html>\n<head>\n<style>\ntable, th, td { border: 1px solid black; border-collapse: collapse; padding: 2px; }\n</style>\n</head>\n"
  stab += "<body>\n"
  stab += "<h1>Table V</h1>\n"
  stab += "<table>\n<tr>\n";
  for n in range(4):
    stab += "<th><i>D<sub>n</sub></i></th>\n"
    stab += "<th><i>B<sub>n</sub></i>/<i>B</i><sub>2</sub><sup><i>n</i>&minus;-1</sup></th>\n"
  stab += "</tr>\n\n"
  stab += mkrow(8,  13, 16)
  stab += mkrow(9,  13, 16)
  stab += mkrow(9,  17, 20)
  stab += mkrow(10, 13, 16)
  stab += mkrow(10, 17, 20)
  stab += "</table>\n"
  stab += "</body>\n"
  stab += "</html>\n"
  return stab



if __name__ == "__main__":
  s = mkTableV()
  open("TableV.html", "w").write(s)
