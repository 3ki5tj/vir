#!/usr/bin/env python



import glob, os, sys, re
import scifmt



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



def getnum(x, err, errmax = 30):
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



def getpct(x, err):
  s = scifmt.scifmt(100*x, err).html(errmax = 10)
  i = s.find("(")
  j = s.find(")")
  # the first character is '+'
  return s[1:i] + s[j+1:]



def mkrow(dim):
  ''' write an HTML row '''

  srow = "<tr>\n<td>%s</td>\n" % dim

  # exact result
  fn = "../../gaussf/gaussfD%smpf.dat" % dim
  Bx, errx = loadvir(fn)

  # DSC
  a = glob.glob("*GFBnPYcD%dn*.dat" % dim)
  if a:
    if len(a[0]) < len(a[1]):
      fn1, fn2 = a[0], a[1]
    else:
      fn1, fn2 = a[1], a[0]
    B1, err1 = loadvir(fn1)
    B2, err2 = loadvir(fn2)
    m = re.search(r"c([0-9\.]+)L(.*)ldbl", fn2)
    kappa = float( m.group(1) )
    lval = int( m.group(2) )

    # compute the error
    em1 = em2 = 0
    for n in range(4, 13):
      em1 = max(em1, abs(B1[n]/Bx[n] - 1))
      em2 = max(em2, abs(B2[n]/Bx[n] - 1))

  for n in range(9, 13):
    srow += "<td style='text-align: left;'>"

    if a:
      # DSC, kappa = 0
      srow += getnum(B1[n], abs(B1[n] - Bx[n])) + "<br>"

      # DSC, kappa != 0
      srow += getnum(B2[n], abs(B2[n] - Bx[n])) + "<br>"

    # exact result
    srow += getnum(Bx[n], errx[n])

    srow += "</td>\n"

  ae = [ 0,   0,    0,    0,     0,     0,
            0.1, 0.01, 0.01, 0.001, 0.001]
  if a: # add a column for kappa
    srow += "<td style='text-align: right;'>"
    srow += "0<br>"
    srow += "(%s)<sub><i>n</i> &ge; %s</sub><br>" % (kappa, lval + 2)
    srow += "&nbsp;"
    srow += "</td>\n"
    srow += "<td style='text-align: right;'>"
    srow += "%s%%<br>%s%%<br>" % ( getpct(em1, ae[dim]), getpct(em2, ae[dim]) )
    srow += "(exact)</td>\n"
  else:
    srow += "<td>&nbsp;</td>\n"
    srow += "<td style='text-align: right;'>(exact)</td>\n"

  srow += "</tr>\n\n"
  return srow



def mkTableIV():
  stab = "<html>\n<head>\n<style>\ntable, th, td { border: 1px solid black; border-collapse: collapse; padding: 2px; }\n</style>\n</head>\n"
  stab += "<body>\n"
  stab += "<h1>Table IV</h1>\n"
  stab += "<table>\n";
  stab += "<tr>\n<th><i>D</i></th>\n"
  for n in range(9, 13):
    stab += "<th><i>B</i><sub>%d</sub>/<i>B</i><sub>2</sub><sup>%d</sup></th>\n" % (n, n - 1)
  stab += "<th><i>&kappa;<sub>n</sub></i></th>\n"
  stab += "<th>Error</th>\n"
  stab += "</tr>\n\n"
  for n in range(1, 11):
    stab += mkrow(n)
  stab += "</table>\n"
  stab += "</body>\n"
  stab += "</html>\n"
  return stab



if __name__ == "__main__":
  s = mkTableIV()
  open("TableIV.html", "w").write(s)
