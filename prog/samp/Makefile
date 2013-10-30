# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm

# optimized compiling
# if -xHost -O3 is used, also define -DRVN_ALIGN=2 or 4
# otherwise it causes trouble for mcgcr2.c in D = 15 and odd dimensions
# for some icc versions
CX = icc -O2 -march=native -ipo \
     -Wall -Wremarks -Wcheck -wd11001 -wd11006
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX =

# OpenMP version
MPCC = $(CX) -openmp
MPLM =

# MPI compiling
MPICC = mpicc -DMPI -Wall -Wextra -O3 -march=native
MPILM =


deps = dg.h dgmapl.h dgcsep.h dgrjw.h dgsc.h dgring.h dgaut.h dghash.h zcom.h Makefile mcutil.h
bins = mcrat0 mcrat1 mcgc2 mcgcr2
bins_d = $(patsubst %, %_d, $(bins))
bins_mp = $(patsubst %, %_mp, $(bins))
bins_mpi = $(patsubst %, %_mpi, $(bins))

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

sync::
	rsync -avL *.[ch] cz1@129.109.88.204:lwork/vir
	rsync -avL *.py cz1@129.109.88.204:lwork/vir
	rsync -avL *.ma cz1@129.109.88.204:lwork/vir
	rsync -avL *.[ch] oo1@stampede.tacc.utexas.edu:work/vir
	rsync -avL *.py oo1@stampede.tacc.utexas.edu:work/vir
	rsync -avL Bring.dat oo1@stampede.tacc.utexas.edu:work/vir

# files are kraken are automatically on keeneland
sync2::
	rsync -e gsissh -avL *.[ch] kraken-gsi.nics.xsede.org:work/vir
	rsync -e gsissh -avL *.py kraken-gsi.nics.xsede.org:work/vir
	rsync -e gsissh -avL Bring.dat kraken-gsi.nics.xsede.org:work/vir

clean:
	rm -rf $(bins) $(bins_d) $(bins_mp) $(bins_mpi) *~ a.out gmon.out \
	  *.his MTSEED bin obj r[0-9][0-9][0-9][a-z]* mr* Zr*

usb:
	$(MAKE) -C ../.. $@

pack:
	$(MAKE) -C ../.. $@

zip:
	$(MAKE) -C ../.. $@

