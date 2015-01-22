#!/bin/bash

gcc -O3 -g -Wall -DLDBL ybgvir.c -lfftw3l -lm -o ybgvir_l
gcc -O3 -g -Wall -DLDBL kirkvir.c -lfftw3l -lm -o kirkvir_l
gcc -O3 -g -Wall -DLDBL -DFFTW ievir.c -lfftw3l -lm -o ievirfftw_l

export Mpt=8388608

# YBG
./ybgvir_l    -n 14 -M 16777216

# Kirkwoord
./kirkvir_l   -n 14 -M 1048576

# PY
./ievirfftw_l -n 14 -M 16777216

# HNC
./ievirfftw_l -n 14 -M $Mpt --hnc

# Hurst
./ievirfftw_l -n 14 -M $Mpt --hm 0.41718

# Rowlinson 1
./ievirfftw_l -n 14 -M $Mpt --rphi 0.16564

# Rowlinson 2
./ievirfftw_l -n 14 -M $Mpt --iphi 0.16564

# HC
./ievirfftw_l -n 14 -M $Mpt --hcs 0.83436

# MS
./ievirfftw_l -n 14 -M $Mpt --ms

# BPGG
./ievirfftw_l -n 14 -M $Mpt --bpggs 1.83436

# Verlet
./ievirfftw_l -n 14 -M $Mpt --va 0.83436 --vb 1.54751

# MP
./ievirfftw_l -n 14 -M $Mpt -q 0.16564

# RY
./ievirfftw_l -n 14 -M $Mpt --rya 0.17015

# Square
./ievirfftw_l -n 14 -M $Mpt --sqrs 0.16564

