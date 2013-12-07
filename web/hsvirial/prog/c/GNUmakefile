# debugging
CC = gcc -O0 -g -Wall -Wunused -Wextra -Wno-unknown-pragmas
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
MPICC = mpicc -DMPI -Wall -Wextra -O3 -march=native -Wno-unknown-pragmas
MPILM =


deps = dg.h dgmap.h dgmapl.h dgutil.h dgcsep.h dgrjw.h dgsc.h \
       dgring.h dgaut.h dgdb.h dghash.h \
       nau0s.h zcom.h GNUmakefile mcutil.h
bins = mcrat0 mcgc2 mcgcr2
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

clean:
	rm -rf $(bins) $(bins_d) $(bins_mp) $(bins_mpi) *~ a.out gmon.out \
	  *.his MTSEED bin obj r[0-9][0-9][0-9][a-z]* mr* Zr*

vircode.zip: $(deps) Bring.dat Z.dat
	zip vircode.zip *.[ch] *.ma *.dat *.py README.usr GNUmakefile

zip: vircode.zip

