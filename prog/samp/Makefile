# debugging
CC = gcc -O0 -g -Wall -Wunused -Wextra -Wno-unknown-pragmas
LM = -lm

# optimized compiling
# if -xHost -O3 is used, also define -DRVN_ALIGN=2 or 4
# otherwise it causes trouble for mcgcr2.c in D = 15 and odd dimensions
# for some icc versions
CX = icc -O2 -ipo \
     -Wall -Wremarks -Wcheck -wd11001 -wd11006 -wd3280
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX =

# OpenMP version
MPCC = $(CX) -openmp
MPLM =

# MPI compiling
MPICC = mpicc -DMPI -Wall -Wextra -O3 -march=native -Wno-unknown-pragmas
MPILM =


deps = dg.h dgmap.h dgmapl.h dgutil.h dgcsep.h dgrjw.h dgsc.h \
       dgring.h dgaut.h dgdb.h dghash.h dgcryr.h \
       nau0s.h zcom.h Makefile mcutil.h
bins = mcrat0 mcgc2 mcgcr2 dbtool
bins_d = $(patsubst %, %_d, $(bins))
bins_mp = $(patsubst %, %_mp, $(bins))
bins_mpi = $(patsubst %, %_mpi, $(bins))
xbins = mds draw

exe: $(bins) $(xbins)

all: $(bins) $(bins_d) $(bins_mp) $(bins_mpi)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bins_mp) : %_mp: %.c $(deps)
	$(MPCC) -o $@ $< $(MPLM)

$(bins_mpi) : %_mpi: %.c $(deps)
	$(MPICC) -o $@ $< $(MPLM)

mds: mds.c
	$(CX) -o $@ $< $(LMX)

draw: draw.c
	$(CC) -DIMLIB2=1 -o $@ $< $(LM) -lglut -lImlib2

%: %.c
	$(CC) -o $@ $< $(LM)

he = cz1@129.109.88.204

synche:
	rsync -avL *.[ch] $(he):lwork/vir
	rsync -avL *.py $(he):lwork/vir
	rsync -avL *.ma $(he):lwork/vir
	rsync -avL Z.dat $(he):lwork/vir

stamp = oo1@stampede.tacc.utexas.edu

sync:
	rsync -avL *.[ch] $(stamp):/work/02464/oo1/vir
	rsync -avL *.py $(stamp):/work/02464/oo1/vir
	rsync -avL Bring.dat $(stamp):/work/02464/oo1/vir
	rsync -avL Z.dat $(stamp):/work/02464/oo1/vir

smic = oo1@smic.hpc.lsu.edu

syncsmic:
	rsync -avL *.[ch] $(smic):work/vir
	rsync -avL *.py $(smic):work/vir
	rsync -avL Bring.dat $(smic):work/vir
	rsync -avL Z.dat $(smic):work/vir

# files are kraken are automatically on keeneland
kraken = kraken-gsi.nics.xsede.org

sync2:
	rsync -e gsissh -avL *.[ch] $(kraken):work/vir
	rsync -e gsissh -avL *.py $(kraken):work/vir
	rsync -e gsissh -avL Bring.dat $(kraken):work/vir
	rsync -e gsissh -avL Z.dat $(kraken):work/vir

clean:
	rm -rf $(bins) $(bins_d) $(bins_mp) $(bins_mpi) $(xbins) *~ a.out gmon.out \
	  mds*.pos mds*.png D*n*.pos \
	  *.his MTSEED bin obj r[0-9][0-9][0-9][a-z]* mr* Zr*

usb:
	$(MAKE) -C ../.. $@

pack:
	$(MAKE) -C ../.. $@

zip:
	$(MAKE) -C ../.. $@

