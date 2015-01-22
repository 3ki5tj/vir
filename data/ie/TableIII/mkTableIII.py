#!/usr/bin/env python



import glob, os, sys
import scifmt



def loadvir(fn, B = None, err = None):
  if not B:   B   = [0]*132
  if not err: err = [0]*132
  for s in open(fn).readlines():
    if s.startswith("#"): continue
    a = s.split()
    n = int(a[0])
    Bn = float(a[1])
    errn = float(a[2])
    if B[n] == 0 or err[n] > errn:
      B[n] = Bn
      err[n] = errn
  return B, err



def getnum(x, err, errmax = 50):
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
      if it > 0:
        s = s[:j] + s[jj+1:]
      break
    else:
      err *= 10
      errmax = 10
    it += 1
  return s



def mkrow(dim, errmax = 10, errmax2 = 50):
  ''' write an HTML row '''

  srow = "<tr>\n<td>%s</td>\n" % dim

  # integral equation
  fn = glob.glob("../xBnPYcD%sn*.dat" % dim)[0]
  B, err = loadvir(fn)

  # Mayer sampling
  fn2ls = glob.glob("../../BnD%sn*.dat" % dim)
  B2, err2 = loadvir(fn2ls[0])
  for fn2 in fn2ls[1:]:
    B2, err2 = loadvir(fn2, B2, err2)

  for n in range(32, 129, 32):
    srow += "<td style='text-align: left;'>"

    # integral equation result
    srow += getnum(B[n], err[n], errmax) + "<br>"

    # Mayer sampling result
    srow += getnum(B2[n], err2[n], errmax2)

    srow += "</td>\n"

  srow += "</tr>\n\n"
  return srow



def mkTableIII():
  stab = "<html>\n<head>\n<style>\ntable, th, td { border: 1px solid black; border-collapse: collapse; }\n</style>\n</head>\n"
  stab += "<body>\n"
  stab += "<h1>Table III</h1>\n"
  stab += "<table>\n";
  stab += "<tr>\n<th><i>D</i></th>\n"
  for n in range(32, 129, 32):
    stab += "<th><i>B</i><sub>%d</sub>/<i>B</i><sub>2</sub><sup>%d</sup></th>\n" % (n, n - 1)
  stab += "</tr>\n\n"
  errmax = 10
  errmax2 = 50
  for n in range(10, 31):
    stab += mkrow(n, errmax, errmax2)
  stab += "</table>\n"
  stab += "</body>\n"
  stab += "</html>\n"
  return stab



if __name__ == "__main__":
  s = mkTableIII()
  open("TableIII.html", "w").write(s)
