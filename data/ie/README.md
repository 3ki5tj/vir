Directories
==============

  * iebenchmark   benchmark data set
  * gaussf        integral equations applied to the Gaussian model
  * pyhc          a set of data with the cavity route computed (mainly for HNC, Fig. 1)
  * lamc          rho-dependent lambda
  * kappa         effect of the shift, applied on hard spheres
  * LJ            Lennard-Jones data set
  * thermo        unsorted thermodynamic data



D               | Method    | Number of points    | Precision*
----------------|-----------|---------------------|-------------------------
2               | DHT       |   524288            |   53-bit (double)
                |           |   262144            |   64-bit (long double)
8,10,...,30     | DHT       |   131072            |   64-bit (long double)
3,5,...,15      | FFT       |  4194304            |   113-bit (quadruple)
17,19,...,27    | FFT       |   262144            |   256-bit
29              | FFT       |   262144            |   384-bit

* Precision is the number of bits to represent the significant.
  The 64-bit double-precision floating-point number has 53 bits.
  The 80-bit extended-precision floating-point number has 64 bits.
  The 128-bit quadruple-precision floating-point number has 113 bits.


Notes
======

D = 20, DHT
  with -M 32768
    double / long double makes little difference
    the final difference at B128 is around 0.007%.
    The multiplier fcorr however differs by 0.25%.
  with -M 65536  the result is strange.
    the multiplier 'fcorr' diverges in the current data file,
    maybe a precision problem?

D = 28, DHT
  with -M 65536 the result diverges at B70-B73.
  to verify the precision.
  should be long double.

D = 26, DHT, nmax = 128
  should use long double, for |Bc| ~ 1e-300 after B80

D = 24, DHT, nmax = 128
  should use long double, for |Bc| ~ 1e-300 after B103

D = 22, DHT -n 128
  fails at B128
