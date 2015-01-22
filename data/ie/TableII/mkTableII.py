#!/usr/bin/env python



import glob, os, sys
import scifmt



def loadvir(fn):
  Bc = [0]*15
  Bv = [0]*15
  By = [0]*15
  for s in open(fn).readlines():
    if s.startswith("#"): continue
    a = s.split()
    n = int(a[0])
    Bc[n] = float(a[1])
    Bv[n] = float(a[2])
    if fn.startswith("x"):
      By[n] = float(a[len(a)/2])
    else:
      By[n - 1] = float(a[-1])
  return Bc, Bv, By



def getnum(x, err):
  while 1:
    s = scifmt.scifmt(x, err).html(errmax = 10, xpmin = -2)
    i = s.find("(")
    j = s.find(")")
    s = s[:i] + s[j+1:]
    k = s.find("&times;")
    if k >= 0:
      s1 = s[:k]
    else:
      s1 = s
    if s1.startswith("&minus;"):
      s1 = s1[6:]
    if len(s1) <= 5:
      err /= 10.0
    else:
      break
  if s.startswith("+"): s = s[1:]
  return s



def mkrow(name, tag, doy):
  ''' write an HTML row '''

  srow = "<tr>\n<td>" + name + "</td>\n<td><i>c</i><br><i>v</i>"
  if doy: srow += "<br><i>y</i>"
  srow += "</td>\n"

  fn = glob.glob("xBn" + tag + "D3*.dat")
  if len(fn) > 0:
    fn = fn[0]
  else:
    fn = glob.glob("Bn" + tag + "D3*.dat")[0]
  Bc, Bv, By = loadvir(fn)

  err = [0, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,
            1e-6, 1e-6, 1e-6, 1e-6, 1e-7,
            1e-7, 1e-7, 1e-7, 1e-7]
  for n in range(4, 13):
    srow += "<td style='text-align: right;'>"
    srow += getnum(Bc[n], err[n])
    srow += "<br>" + getnum(Bv[n], err[n])
    if doy:
      srow += "<br>" + getnum(By[n], err[n])
    srow += "</td>\n"

  srow += "</tr>\n\n"
  return srow



def mkTableII():
  stab = "<html>\n<head>\n<style>\ntable, th, td { border: 1px solid black; border-collapse: collapse; padding: 2px 4px 2px 4px; }\n</style>\n</head>\n"
  stab += "<body>\n"
  stab += "<h1>Table II</h1>\n"
  stab += "<table>\n";
  stab += "<tr>\n<th>Integral equation</th>\n<th>&nbsp;</th>\n"
  for n in range(4, 13):
    stab += "<th><i>B</i><sub>%d</sub>/<i>B</i><sub>2</sub><sup>%d</sup></th>\n" % (n, n - 1)
  stab += "</tr>\n\n"
  stab += mkrow("YBG",          "YBG",  True)
  stab += mkrow("Kirkwood",     "K",    True)
  stab += mkrow("PY",           "PY",   True)
  stab += mkrow("HNC",          "HNC",  True)
  stab += mkrow("Hurst",        "H",    False)
  stab += mkrow("Rowlinson 1",  "R",    False)
  stab += mkrow("Rowlinson 2",  "IR",   False)
  stab += mkrow("HC",           "HC",   False)
  stab += mkrow("MS",           "MS",   False)
  stab += mkrow("BPGG",         "BPGG", False)
  stab += mkrow("Verlet",       "V",    False)
  stab += mkrow("MP",           "MP",   False)
  stab += mkrow("RY",           "RY",   False)
  stab += mkrow("Quadratic",    "SQR",  False)
  stab += "</table>\n"
  stab += "</body>\n"
  stab += "</html>\n"
  return stab



if __name__ == "__main__":
  s = mkTableII()
  open("TableII.html", "w").write(s)
