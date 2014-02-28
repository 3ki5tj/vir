# debugging
CC = gcc -O0 -g -Wall -Wunused -Wextra -Wno-unknown-pragmas
LM = -lfftw3 -lm

# optimized compiling
# if -xHost -O3 is used, also define -DRVN_ALIGN=2 or 4
# otherwise it causes trouble for mcgcr2.c in D = 15 and odd dimensions
# for some icc versions
CX = icc -O2 -ipo \
     -Wall -Wremarks -Wcheck -wd11001 -wd11006
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX = -lfftw3

# OpenMP version
MPCC = $(CX) -openmp
MPLM = -lfftw3_omp

# MPI compiling
MPICC = mpicc -DMPI -Wall -Wextra -O3 -march=native -Wno-unknown-pragmas
MPILM = -lfftw3_mpi


csrc = $(wildcard *.c)
deps =
bins = $(patsubst %.c, %, $(csrc))
bins_d = $(patsubst %, %_d, $(bins))

exe: $(bins)

all: $(bins) $(bins_d) $(bins_mp) $(bins_mpi)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bins_mp) : %_mp: %.c $(deps)
	$(MPCC) -o $@ $< $(MPLM)

$(bins_mpi) : %_mpi: %.c $(deps)
	$(MPICC) -o $@ $< $(MPLM)

syncipo:
	rsync -avL *.[ch] cz1@129.109.88.204:lwork/vir

sync: syncipo
	rsync -avL *.[ch] oo1@stampede.tacc.utexas.edu:/work/02464/oo1/vir

# files are kraken are automatically on keeneland
sync2::
	rsync -e gsissh -avL *.[ch] kraken-gsi.nics.xsede.org:work/vir

clean:
	rm -rf $(bins) $(bins_d) $(bins_mp) $(bins_mpi) *~ a.out gmon.out \
	  *.his MTSEED bin obj r[0-9][0-9][0-9][a-z]* mr* Zr*

usb:
	$(MAKE) -C ../.. $@

pack:
	$(MAKE) -C ../.. $@

zip:
	$(MAKE) -C ../.. $@

