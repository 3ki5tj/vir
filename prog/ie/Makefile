# debugging
CC = gcc -O0 -g -Wall -Wunused -Wextra -Wno-unknown-pragmas
LM = -lm

# optimized compiling
# warnings 11001 and 11006 are notices invoked by -ipo
CX = icc -O2 -ipo -Wall -Wremarks -Wcheck -wd11001 -wd11006
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX =

LFFTW = -lfftw3
LFFTWQ = -lfftw3q -lquadmath
LMPFR = -lmpfr -lgmp
LGMP = -lgmp
LGSL = -lgsl -lgslcblas

# source code that uses GMP and MPFR
csrc            = $(wildcard *.c)
csrc_fftw       = $(wildcard *fftw*.c)
csrc_mpfr       = $(wildcard *mpfr*.c)
csrc_gmp        = $(wildcard *gmp*.c)
csrc_gsl        = $(wildcard *gsl*.c)
csrc_extra      = $(filter-out \
		    $(csrc_fftw) $(csrc_mpfr) $(csrc_gmp) $(csrc_gsl), \
		    $(csrc))

deps = Makefile ieutil.h ieutilmpfr.h fft.h mpfft.h
bins            = $(patsubst %.c, %,   $(csrc))
bins_q          = $(patsubst %.c, %_q, $(csrc))
bins_d          = $(patsubst %,   %_d, $(bins))
bins_fftw       = $(patsubst %.c, %,   $(csrc_fftw))
bins_fftw_q     = $(patsubst %.c, %_q, $(csrc_fftw))
bins_fftw_d     = $(patsubst %,   %_d, $(bins_fftw))
bins_mpfr       = $(patsubst %.c, %,   $(csrc_mpfr))
bins_mpfr_d     = $(patsubst %,   %_d, $(bins_mpfr))
bins_gmp        = $(patsubst %.c, %,   $(csrc_gmp))
bins_gmp_d      = $(patsubst %,   %_d, $(bins_gmp))
bins_gsl        = $(patsubst %.c, %,   $(csrc_gsl))
bins_gsl_d      = $(patsubst %,   %_d, $(bins_gsl))
bins_extra      = $(patsubst %.c, %,   $(csrc_extra))
bins_extra_d    = $(patsubst %,   %_d, $(bins_extra))

exe: $(bins)

all: $(bins) $(bins_d) $(bins_q)

$(bins_fftw) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LFFTW)

$(bins_fftw_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LFFTW)

$(bins_fftw_q) : %_q: %.c $(deps)
	$(CC) -DQUAD -o $@ $< $(LM) $(LFFTWQ)

$(bins_mpfr) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LMPFR)

$(bins_mpfr_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LMPFR)

$(bins_gmp) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LGMP)

$(bins_gmp_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LGMP)

$(bins_gsl) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LGSL)

$(bins_gsl_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LGSL)

$(bins_extra) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_extra_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

syncipo:
	rsync -avL *.[ch] cz1@129.109.88.204:lwork/vir

sync: syncipo
	rsync -avL *.[ch] oo1@stampede.tacc.utexas.edu:/work/02464/oo1/vir

# files are kraken are automatically on keeneland
sync2::
	rsync -e gsissh -avL *.[ch] kraken-gsi.nics.xsede.org:work/vir

clean:
	rm -rf $(bins) $(bins_d) $(bins_q) \
	  *~ a.out gmon.out MTSEED bin obj \
	  r[0-9][0-9][0-9][a-z]* vv*.dat

usb:
	$(MAKE) -C ../.. $@

pack:
	$(MAKE) -C ../.. $@

zip:
	$(MAKE) -C ../.. $@
