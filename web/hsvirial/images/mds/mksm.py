import glob, os

fns = glob.glob("*.png")

for fn in fns:
  a, b = os.path.splitext(fn)
  fnsm = a + "sm" + b
  os.system("convert -resize %sx%s %s %s" % (
    240, 240, fn, fnsm))
