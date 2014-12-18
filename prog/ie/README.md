# Compute virial coefficients from integral equations #

## List of files ##

### Core library ###

  Code        |  Dependencies | Description
------------- | ------------- | --------------
  zcom.h      |   none        | common routines
  xdouble.h   |   none        | defines the datatype `xdouble`, which can be `double`, `long double` or `__float128`.
  ieutil.h    |   none        | common routines for integral equations
  fft.h       |   none        | home-made code for FFT, supports `xdouble.h`, used as a reference implementation
  slowdht.h   |   GSL         | discrete Hankel transform (DSC), supports `xdouble.h`), permits disk-based Hankel transformation to support transforms of many evaluation points
  fftx.h      |   FFTW3/GSL   | a wrapper of fft.h and slowdht.h


### For virial coefficients ###

#### C programs (native datatypes family) ####

  Code        | Dependencies  | Description
------------- | ------------- | --------------
  ievir.c     | FFTW3/GSL     | compute virial coefficients of integral equations
  iegsl.c     |   GSL         | alias of `ievir.c` as a symblic link, for even/odd dimensions
  ieodfftw.c  |  FFTW3        | alias of `ievir.c` as a symbolic link, for odd dimensions
  snapfixgsl.c (deprecated) | GSL | fix the snapshotXXX_.dat output of iegsl.c
  ozcrtr.c    | FFTW3/GSL     | given c(r), compute t(r) from the OZ relation
  ljievir.c   | GSL/FFTW3     | for the LJ potential (special PYA closure)


#### C programs (MPFR family) ####

  Code           |  Description
---------------- | ----------------------
  ieodmpfr.c     | generalized odd-dimensional case
  ieutilmpfr.h   | common routines for the MPFR programs, the counterpart of ieutil.h
  fftmpfr.h      | multiple-precision code for FFT, the counterpart of fft.h

#### Python script ####

  Code           | Description
---------------- | -----------------------
  ievirextra.py  | extrapolate virial coefficients for different grid spacings
  ieoptparam.py  | script to find optimal parameters of kappa in the detailed self-consistent condition
  scifmt.py      | format a number with error in parentheses used in ievirextra.py

#### Additional code ####

See bak/README for additional code.
ie3dxxx.c is kept along with ieodxxx.c for the purpose of double checking.



### Integral equation (expansion around a density) ###

  Code          | Dependencies  | Description
--------------- | ------------- | -------------
  ie.c          | FFTW3/GSL     | integral equation solver
  ielmv.h       | none          | LMV solver
  iemdiis.h     | none          | MDIIS solver





## Programs usage ##

### ievir.c ###

#### Overview ####

  * compute virial coefficients from a specific closure

  * works for any dimensions
    - For even dimensions, define -DDHT to use the discrete Hankel transform (DHT) during compiling.
      This version has the alias `iegsl.c`, which is a symbolic link to `ievir.c`
      because it uses the GNU Scientific Library (GSL for the Bessel functions.)

    - For odd dimensions, define -DFFTW to use FFTW3 during compilation.
      This version has the alias `ieodfftw.c`, which is symbolic link to `ievir.c`

  * supports the following closures
    - Percus-Yevick (PY)
    - hypernetted-chain (HNC)
    - Rowlinson, t(r) = (y(r) - 1) (1 - phi) + phi ln y(r)
    - inverse Rowlinson (IR), y(r) = exp(t(r)) phi + (1 - phi) (1 + t(r))
    - Hurst, t(r) = y(r)^m log y(r)
    - Hutchinson-Conkie (HC)
    - Verlet (VM)
    - Rogers-Young (RY)
    - Marucho-Pettitt (MP)
    - Detailed (term-by-term) self-consistent (DSC)
      This is based on the inverse-Rowlinson closure.
    - Lambda-Detailed self-consistent (&lambda;-DSC)
      This closure tunes the parameter, if any, of the above closures
      to the reach the self-consistency for each order.

  * supports check points, (with option "-s")
    which allows a restartable run.



#### Source code ####

  ievir.c

#### Compilation through the Makefile ####

Enter
```
make iegsl ieodfftw
```
Here,
`iegsl` is mainly for even dimensions,
`ieodfftw` is for odd dimensions.


#### Manual compilation ####

Normal example, for even/odd dimensions
```
gcc -O3 -DDHT ievir.c -lgsl -lgslcblas -lm
```
To run it
```
./a.out -D 20 --corr -n 64 -R 67 -M 32768
```

For odd dimensions only
```
gcc -O3 -DFFTW ievir.c -lfftw3 -lm
```

To `long double`, define `LDBL` or `LONG` in compilation
```
gcc -O3 -DLDBL -DDHT iegsl.c -lgsl -lgslcblas -lm
```

For slowdht with `__float128`, define `F128` or `QUAD` in compilation
```
gcc -O3 -DF128 -DDHT -Wall -Wextra iegsl.c -lgsl -lgslcblas -lquadmath -lm
```
Or for odd dimensions,
```
gcc -O3 -DF128 ieodfftw.c -lfftw3q -lquadmath -lm
```
Test it by
```
./a.out -D 15 -n 64 -R 66.01 -M 262144 --corr
```
Remember to use `gcc` instead of `icc` in this case,
because `__float128` is supported by GCC only.


#### Useful command-line options ####

  Options       | Description
--------------- | ----------------
  -h            | list all options
  -D30          | specifiy the dimension as 30
  -n128         | specifiy the maximal order as 128
  -M8192        | specify 8192 evaluation points along r
  --disk=2      | use disk-based discrete Hankel transform
  --hnc         | do the HNC closure instead of the PY one
  --corr        | do the detailed self-consistent (DSC) closure
  --lamc        | do the &lambda;-DSC

### Precision ###

Tested for dimensions 2 - 30.

Generally, we need more precise data type for higher dimensional cases.

The following examples may help the user to choose the proper parameters

#### Example 1 ####

For D = 30, n = 128, M = 3072 (sampling points)

  PY:   __float128 and long double make no difference
        up to 14 significant digits.

        For B(128):
                              Bc                    Bv
          long double:     -1.65055622495195e7    -1.12865683738071e8
          __float128:      -1.65055622495195e7    -1.12865683738071e8
        So there is a 20% difference in Bc.

        By comparison, with -M2048, B(128) are given below
                              Bc                    Bv
          long double:     -3.14733368059798e1    -2.34259117610089e3
          __float128:      -3.14733368059798e1    -2.34259117610089e3

        However, the __float128 version is 12 times
        as slow as the long double version.
        (1400s vs 115s, disk-based, no matrix time, T60).

  HNC:  The __float128 and long double versions make
        no difference in Bv, but significant difference in Bc
        For B(128):
                              Bc          Bv
          long double:     -8.0428e7    -1.11432807771606e8
          __float128:      -1.0507e7    -1.11432807771606e8
        So there is a 20% difference in Bc.

        By comparison, with -M2048, B(128) are given below
                              Bc          Bv
          long double:     -1.508e3     -2.31117866186384e3
          __float128:      -2.004e3     -2.31117866186384e3

        The __float128 version is much slower than
        the long double version.
        (1502 vs 61s, disk-based, no matrix time, T60)

  SC:   SC stands for the Self-Consistent closure.
        The __float128 version gives more accurate results
        than long double version.
        For B(128):
          long double:    -1.1138e+8
          __float128:     -1.1154e+8
        So the difference is 0.142%.

        By comparison, with -M2048, B(128) are given below
          long double:    -2.3072e+3
          __float128:     -2.3105e+3
        So the difference is also about 0.157%.
        It is possible that the difference diminishes with
        the number of evaluation points.

        Note that the figure themselves are very wrong,
          the correct figure from -M65536 (long double) is -6.8e11.
        Thus the bin-size dependence is much larger than
          the dependence on the float-point precision.

        The __float128 version is much slower than
        the long double version.
        (1497 vs 61s, disk-based, no matrix time, T60)

        Precision of using long double and double
        -M4096 -n128 compared with __float128
        see data/iebenchmark
        --------------------------------------
            D       long double       double
            2         4.76e-7         2.53e-5
                      8.88e-16        2.75e-10  (-n32)
            4         1.04e-4         7.36e-2
                      2.97e-11        5.15e-14  (-n16)
            6         2.22e-16        9.73e-14
                      0               1.58e-14  (-n64)
            8         0               1.09e-13
            10        3.55e-15        3.62e-12
            12        8.04e-14        8.31e-11
            14        2.28e-13        1.15e-8
            16        1.81e-10        2.87e-7
            18        4.02e-9         2.46e-6
            20        7.64e-8         5.71e-5
            22        1.29e-6         1.75e-3
            24        1.97e-5         underflow
            26        7.53e-5         underflow
            28        5.30e-4         underflow
            30        1.26e-3         underflow



## ieodmpfr ##


Overview
---------
High precision version of ieodfftw.c (ievir.c), currently with limited features.


Compilation
------------

When `ieodfftw` (`ievir.c` compiled with `-DFFTW` fails,
switch to the mpfr version to enable higher precisions.
```
icc ieodmpfr.c -lmpfr -lgmp && ./a.out -D 15 -n 128 -R 131.072 -M 262144 --corr
```
Note, FFTW is unavailable in this case, and we use `fftmpfr.h` in this case.

