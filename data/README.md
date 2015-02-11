# Data directory #

## Files ##

File        | Description
------------|-------------------------------------------------
pack.sh     |  pack data for Bossman
Bring.dat   |  values of the Ring integral computed from the Mathematica script `prog/ring.ma`
Z.dat       |  compiled values of the partition function
mkZ.py      |  script to compile `Z.dat`
updateZ.sh  |  shell script that collects data run on kraken, then assembles them by `mkZ.py`
mkeos.py    |  equation of state script (unfinished)
fbtime.dat  |  time to compute fb



Directory   | Description
------------|-------------------------------------------------
pyhs        |  analytic PY solution (`prog/pyhs`)
pydata      |  a duplicate of `pyhs`, for making copies to Bossman
ie          |  integral equation data (`prog/ie`)
gaussf      |  Gaussian model data (`prog/gaussf`)
stampede    |  data collected on `stampede.tacc.utexas.edu`
kraken      |  data collected on `kraken-gsi.nics.xsede.org`
supermic    |  data collected on `smic.hpc.lsu.edu`
LJ          |  standard Lennard-Jones data
