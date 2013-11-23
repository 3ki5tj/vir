#!/usr/bin/env

''' make Bring.dat Java array '''

lines = open("Bring.dat").readlines()

src = "  double [][] BringArr = new double [][] {\n"
src += "    {" + ", ".join(["1",]*64) + "},\n"
src += "    {" + ", ".join(["1",]*64) + "},\n"
for line in lines:
  line = line.strip()
  if line == "": continue
  arr = line.split()
  dim = int(arr[0])
  arr[0] = "1";
  src += "    {" + ", ".join(arr) + "},\n"
src = src.rstrip(",\n")
src += "};\n"
print src
