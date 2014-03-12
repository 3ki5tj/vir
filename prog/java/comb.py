#!/usr/bin/env python



''' Combine Java source code '''

proj = "VirSamp"

ls = [proj + "App.0.java",
      "MCSamp.java",
      "Diagram.java",
      "Bits.java",
      "Ave.java",
      "MyScheme.java",
      "MyCanvas.java",
      "XYZModel.java",
      "Atom.java",
      "Matrix3D.java",
      "MDS.java"]

def trim(s):
  for i in range(len(s)):
    if s[i].startswith("import "):
      s[i] = ""

src = []
for fn in ls:
  s = open(fn).readlines()
  if not fn.startswith(proj): trim(s)
  src += s

fnout = proj + "App.java"
open(fnout, "w").writelines(src)

