#!/bin/bash

gcc -O3 -g -Wall        -DFFTW ie.c -lfftw3  -lm -o iefftw
gcc -O3 -g -Wall -DLDBL -DFFTW ie.c -lfftw3l -lm -o iefftw_l


# SMSA, compressibility route
./iefftw_l --LJ=split -Smdiis -R163.84 -N 32768 -T 1.071 --rho=0.3107 --cr=smsa_c.dat -v

# SMSA, virial route
./iefftw_l --LJ=split -Smdiis -R163.84 -N 32768 -T 1.491 --rho=0.30555 --cr=smsa_v.dat -v

# PY, compressibility route
./iefftw_l -Smdiis -R163.84 -N 32768 -T 1.3199 --rho=0.278 --cr=py_c.dat -v

# PY, virial route (not precise)
./iefftw_l -Smdiis -R163.84 -N 32768 -T 1.32 --rho=0.2503 --cr=py_v.dat -v

# self-consistent
./iefftw --LJ=split -R163.84 -N65536 --LJlrs=0.3061 -T 1.2998 --rho=0.2970 --cr=sc_c.dat -v

