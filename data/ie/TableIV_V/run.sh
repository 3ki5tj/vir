#!/bin/bash

gcc -O3 -g -Wall -DLDBL -DDHT   ievir.c -lgsl -lgslcblas            -lm -o ievirgsl_l
gcc -O3 -g -Wall -DLDBL -DFFTW  ievir.c -lfftw3l                    -lm -o ievirfftw_l
#gcc -O3 -g -Wall -DQUAD -DDHT   ievir.c -lgsl -lgslcblas -lquadmath -lm -o ievirgsl_q
#gcc -O3 -g -Wall -DQUAD -DFFTW  ievir.c -lquadmath -lfftw3q         -lm -o ievirfftw_q

./ievirgsl_l -G -D6 -n 20 -M 8192 --corr
./ievirgsl_l -G -D6 -n 20 -M 8192 --corr -c 0.047 -L6

./ievirfftw_l -G -D7 -n 14 -M 65536 --corr
./ievirfftw_l -G -D7 -n 14 -M 65536 --corr -c 0.0085 -L6

./ievirgsl_l -G -D8 -n 20 -M 8192 --corr
./ievirgsl_l -G -D8 -n 20 -M 8192 --corr -c 0.00325 -L6

./ievirfftw_l -G -D9 -n 24 -M 16384 --corr
./ievirfftw_l -G -D9 -n 24 -M 16384 --corr -c 0.000924 -L5

./ievirgsl_l -G -D10 -n 24 -M 8192 --corr
./ievirgsl_l -G -D10 -n 24 -M 8192 --corr -c 0.00029 -L4

