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
scifmt.py   |  print a number with error in scientific format, e.g. 3.1415(21)
virsum.py   |  used to produce summary file for each dimension, BnDxnyyy.dat


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


## Ring Integrals, Bring.dat
Ring integrals are the reference integrals for computing the virial coefficients.
Ring integral can be analytically computed.

## Partition Functions, Z.dat
Partition functions are another set of reference integrals for computing the virial coefficients.
The partition is always positive, and it encompasses are cluster configurations relevant to the computation of virial coefficients.  Hence, it serves a better reference than the ring integral.

However, partition functions cannot be analytically computed except in a few simple cases.
It can be computed from grand canonical Monte Carlo simulation, in which the order n, or the number of particles varies.

## Summary Data Files, BnDxnyyy.dat

These are created by the summarizing script virsum.py.
Most of the files are under the stampede directory

BnD11n22.dat contains results for the simulation on the D = 11 dimensional hard-sphere fluid, with the maximal order n = 22.
This is text file, so you can open it by Vim, nano or any Linux/Unix text editor.
* The first column is the order.
* The second column is the virial coefficient in the form of Bn/B2^(n-1).
* The third column is the estimated standard deviation.
* The last column is the total sample size. Note that the samples come from a Monte Carlo simulations, so they are by no means independent.


Now for the 11th dimension, we have another file, BnD11n32.dat, which goes up to n = 32.  This file comes from an independent grand-canonical simulation.  However, it has poorer precision for the lower order (as shown in column three).  So, for a lower-order virial coefficient, I would use the numbers from the former file.

## Zrhxxx.dat and Zrxx.dat
These are created by summarizing script
stampede/sum.py
