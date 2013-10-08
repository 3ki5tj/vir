# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm
# optimized compiling
# if -xHost -O3 is used, also define -DRVN_ALIGN=2 or 4
# otherwise it causes trouble for mcgcr2.c in D = 15 and odd dimensions
# for some icc versions
CX = icc -O2 -march=native -ipo \
     -Wall -Wunused -Wremarks -Wcheck -wd11001 -wd11006
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX =

deps = dg.h dgcsep.h dgrjw.h dgsc.h dgring.h zcom.h Makefile mcutil.h
bins = mcrat0 mcrat mcgc2 mcgcr2
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

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
	rm -rf $(bins) $(bins_d) *~ a.out gmon.out *.his MTSEED bin obj

usb:
	$(MAKE) -C ../.. usb
