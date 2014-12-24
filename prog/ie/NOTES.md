gaussf running
===============
n 13, e 24
n 14, e 23
n 15, e 23
n 17, e 24
n 18, e 24
n 19, e 25
n 20, e 25


Shift of the hard sphere model
==============================

New shift
----------

Hard sphere
~~~~~~~~~~~~

D = 2
* best parameters: ./ievirgsl --corr -D2 -c0.194 -d0.0 -L4, (c' 0.194), error 1.54203833119% (n 6, 1.07464109577%)
* best parameters: ./ievirgsl --corr -D2 -c-0.06 -d0.04 -L3, (c' -0.14), error 0.7131301684% (n 11, 0.468641492886%)

D = 3
* best parameters: ./ievirfftw_q -M8192 --corr -D3 -c0.335 -d0.0 -L3, (c' 0.335), error 6.31580626675% (n 11, 3.78196593326%)

D = 6
* best parameters: ./ievirgsl --corr -D6 -c0.352 -d0.0 -L2, (c' 0.352), error 8.28359264602% (n 12, 4.09865752935%)
* best parameters: ./ievirgsl --corr -D6 -c0.45 -d-0.04 -L2, (c' 0.53), error 7.83816390485% (n 5, 4.38635292495%)
fit to n = 16
* best parameters: ./ievirgsl --corr -D6 -c0.36 -d0.0 -L2, (c' 0.36), error 17.5190533199% (n 15, 8.57740113315%)

Gaussian model
~~~~~~~~~~~~~~~

D = 5
* best parameters: ./ievirfftw_q -M8192 --corr -G -D5 -c-0.136 -d0.0 -L3, (c' -0.136), error 12.6855302858% (n 5, 7.83671310393%)

D = 6
* best parameters: ./ievirgsl --corr -G -D6 -c0.047 -d0.0 -L6, (c' 0.047), error 2.58551856198% (n 12, 1.05634911603%)
* best parameters: ./ievirgsl --corr -G -D6 -c-0.04 -d0.011 -L5, (c' -0.062), error 0.479199882519% (n 8, 0.269281246909%)

D = 7
* best parameters: ./ievirfftw_q -M8192 --corr -G -D7 -c0.0085 -d0.0 -L6, (c' 0.0085), error 0.506242267168% (n 12, 0.210889303082%)
* best parameters: ./ievirfftw_q -M8192 --corr -G -D7 -C-0.004 -d0.002 -L3, (c 2.77555756156e-17), error 0.526348985362% (n 6, 0.332090456191%)

D = 8
* best parameters: ./ievirgsl --corr -G -D8 -c0.00325 -d0.0 -L6, (c' 0.00325), error 0.128332699584% (n 7, 0.0565939935881%)

D = 9
* best parameters: ./ievirfftw_q -M8192 --corr -G -D9 -c0.00089 -d0.0 -L5, (c' 0.00089), error 0.0385239318704% (n 12, 0.0209715618515%)




Old shift (before Aug 17, 2014)
------------------------------

2D hard-sphere:
./ievirgsl --corr -D2 -c-0.132 -d0.061 -L3, (c' -0.254), error 0.427821686445% (n 6)
./ievirgsl --corr -D2 -c0.205 -L4, (c' 0.205), error 1.02089359331% (n 10)

3D hard-sphere:
./ievirfftw_q -M8192 --corr -D3 -c-0.374 -d0.204 -L2, (c' -0.782), error 2.88427972437% (n 12)
./ievirfftw_q -M8192 --corr -D3 -c0.506 -d0.0 -L3, (c' 0.506), error 9.76567947096% (n 11)

4D hard-sphere
./ievirgsl --corr -D4 -c-0.675 -d-0.038 -L2, (c' -0.599) error 54.1577541112% (n 5)
./ievirgsl --corr -D4 -c-0.672 -d-0.037 -L2, error 73.5835955182% (n 9)

5D hard-sphere
./ievirfftw_q -M8192 --corr -D5 -c-0.876 -d-0.017 -L2, (c' -0.842), error 214.437488766% (n 8)
./ievirfftw_q -M8192 --corr -D5 -c-0.874 -d-0.018 -L2, (c' -0.838), error 214.466127363% (n 7)

6D hard-sphere
./ievirgsl --corr -D6 -c0.182 -d0.209 -L2, (c' -0.236), error 6.22214169297% (n 7)
./ievirgsl --corr -D6 -c0.675 -d0.0 -L2, error 11.0704659824% (n 8)

Gaussian model
===============

5D gaussian
./ievirfftw_q -M8192 --corr -G -D5 -c-0.103 -d-0.007 -L3, (c' -0.089), error 12.462324563% (n 5)



## Critical point ##

### SMSA/PYA approximation ###

#### Compressibility route ####

rho in the following table is the density
that achieves the zero d^2(beta P)/d(rho)^2

  T     |  rho    |  d(beta P)/d(rho)   |  d^2(beta P)/d(rho)^2
--------+---------+---------------------+------------------------
  1.28  | 0.2645  |  0.305              |  0.01
  1.27  | 0.2662  |  0.290              |  0.00
  1.26  | 0.2682  |  0.276              | -0.001
  1.25  | 0.2703  |  0.2614             |  0.000
  1.20  | 0.281   |  0.1856             |  0.002
  1.15  | 0.292   |  0.1071             |  0.000
  1.10  | 0.304   |  0.0328             |  0.01
  1.09  | 0.306   |  0.0202             | -0.00
  1.08  | 0.3086  |  0.0094             |  0.00
  1.072 | 0.311   |  0.0029             |  0.00
  1.071 | 0.311   |  0.0022             |  0.002

A typical command is the following
```
./iefftw --LJ=split -T 1.08 -S2 -R163.84 -N 32768 --rho=0.305
```

Critical point
T = 1.071, rho = 0.311.


#### Virial route ####

  T     |  rho    |  d(beta P)/d(rho)   |  d^2(beta P)/d(rho)^2
--------+---------+---------------------+------------------------
  1.5   | 0.304   |  0.017              | -0.000
  1.49  | 0.3057  | -0.002              | -0.001
  1.48  | 0.3074  | -0.021              | -0.001
  1.4   | 0.321   | -0.189              | -0.006

Critical point
T = 1.49, rho = 0.3057.


### PY approximation ###

#### Compressibility route ####

  T     |  rho    |  d(beta P)/d(rho)   |  d^2(beta P)/d(rho)^2
--------+---------+---------------------+------------------------
  1.40  | 0.275   |  0.141              |  0.01
  1.35  | 0.280   |  0.055              | -0.01
  1.32  | 0.278   |  0.001              | -0.00

Critical point
T = 1.32, rho = 0.278.


#### Virial route ####

  T     |  rho    |  d(beta P)/d(rho)   |  d^2(beta P)/d(rho)^2
--------+---------+---------------------+------------------------
  1.325 | 0.2456  |  0.106              |  0.00
  1.32  | 0.2503  |  0.09               |  0.00

Critical point
T = 1.xx, rho = 0.3xx


### HNC approximation ###

#### Compressibility route ####

  T     |  rho    |  d(beta P)/d(rho)   |  d^2(beta P)/d(rho)^2
--------+---------+---------------------+------------------------
  1.5   | 0.277   |  0.243              |   -0.002
  1.48  | 0.280   |  0.211              |    0.001
  1.45  | 0.283   |  0.159              |    0.000
  1.42  | 0.283   |  0.097              |   -0.005
        | 0.2831  |  0.097              |   -0.002
        | 0.2832  |  0.097              |   -0.001
  1.41  | 0.280   |  0.066              |    0.029
  1.409 | 0.279   |  0.061              |    0.03
  1.408 | 0.277   |  0.054              |   -0.01
  1.407 | 0.2767  |  0.041              |   -1.29

Critical point
T = 1.4, rho = 0.277


### SC partial PY ###

Critical point
ljlrs = 0.3051, T = 1.2998, rho = 0.296444, bP = 0.083, bP/rho = 0.280

./iefftw --ljlrs=0.3051 --LJ=split -T 1.2998 --rho=0.297 -R81.92 -N65536 --damp=0.1 -v

##### ljlrs = 0.24 #####

T = 1.315, rho = 0.2935,  dbPc/v  0.006/0.05
 
##### ljlrs = 0.25 #####

T = 1.315, rho = 0.2939,  dbPc/v  0.00729/0.04482
 
##### ljlrs = 0.27 #####

T = 1.315, rho = 0.2947,  dbPc/v  0.010/0.0344

##### ljlrs = 0.29 #####

T = 1.315, rho = 0.2953,  dbPc/v  0.013/0.024
 
##### ljlrs = 0.30 #####

T = 1.315, rho = 0.2956,  dbPc/v  0.015/0.018
T = 1.310, rho = 0.29605, dbPc/v  0.009/0.012
T = 1.305, rho = 0.29626, dbPc/v  0.003/0.0066
T = 1.303, rho = 0.29624, dbPc/v  0.001/0.005
 
##### ljlrs = 0.304 #####

T = 1.303, rho = 0.29644, dbPc/v  0.0018/0.0024
 
##### ljlrs = 0.3047 #####

T = 1.303, rho = 0.29647, dbPc/v  0.00190/0.00196
T = 1.302, rho = 0.29645, dbPc/v  0.00118/0.00122
T = 1.301, rho = 0.29643, dbPc/v  0.00058/0.00067
T = 1.300, rho = 0.29643, dbPc/v  0.00010/0.00038
 
##### ljlrs = 0.3051 #####

T = 1.2998, rho = 0.296444, dbPc/v  0.00004/0.00003
 
##### ljlrs = 0.31 #####

T = 1.315, rho = 0.2958,  dbPc/v  0.017/0.013
 
##### ljlrs = 0.32 #####

T = 1.315, rho = 0.29605,  dbPc/v  0.019/0.008
 
##### ljlrs = 0.33 #####

T = 1.315, rho = 0.2962,  dbPc/v  0.021/0.002
 
