#!/usr/bin/env python



''' Combine Java source code '''

proj = "VirSamp"

ls = [proj + "App.0.java",
      "MCSamp.java",
      "Diagram.java",
      "DiagramMap.java",
      "Bits.java",
      "Ave.java",
      "XYScheme.java",
      "XYZCanvas.java",
      "XYZModel.java",
      "XYZModelMC.java",
      "Atom.java",
      "Matrix3D.java",
      "MDS.java"]

def trim(s):
  for i in range(len(s)):
    if s[i].startswith("import "):
      s[i] = ""

src = []
for i in range(len(ls)):
  fn = ls[i]
  s = open(fn).readlines()
  if i > 0: trim(s)
  src += s

fnout = proj + "App.java"
open(fnout, "w").writelines(src)

