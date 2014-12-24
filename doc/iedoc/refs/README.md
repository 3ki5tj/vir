## Zerah1986 paper ##

Table II. the SMSA column can be reproduced
by the program `ie.c` or `iefftw.c`.
Using the first row for example,

./iefftw -T 100 --rho=0.1 --LJ=split


### beta P / rho ###

We could not reproduce the T = 100, rho = 2, 2.5 data points.
The densities are way too large to be physical.

./iefftw -T 20 --rho=1.3333 --LJ=split --damp=0.5
./iefftw -T 20 --rho=1.765 --LJ=split --damp=0.5 -SLMV -M 10

./iefftw -T 5 --rho=0.5 --LJ=split
./iefftw -T 5 --rho=1 --LJ=split
./iefftw -T 5 --rho=1.279 --LJ=split

The data for T = 5, rho = 0.1 should probably be 5.73 instead of 5.27

./iefftw -T 2.74 --rho=0.55 --LJ=split
./iefftw -T 2.74 --rho=1.0 --LJ=split
./iefftw -T 2.74 --rho=1.1 --LJ=split


