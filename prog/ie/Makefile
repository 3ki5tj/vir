prj = ievir

LQ = -lquadmath
LFFTW = -lfftw3
LFFTWL = -lfftw3l
LFFTWQ = -lfftw3q
LMPFR = -lmpfr -lgmp
LGMP = -lgmp
LGSL = -lgsl -lgslcblas

# debugging
CC = gcc -O0 -g3 -Wall -Wunused -Wextra -Wno-unknown-pragmas
LM = -lm

# optimized compiling
# warnings 11001 and 11006 are notices invoked by -ipo
# warning 1572 is the float-point comparison
# use -O1 for compiling speed, use -O2 for speed
CCOPT = -O1
CX = icc $(CCOPT) -ipo -Wall -Wremarks -Wcheck -wd11001 -wd11006 -wd1572 -wd981 -wd3280
#CX = gcc -O3 -march=native -fwhole-program -Wall -Wextra -Wunused
LMX =

# optimized compiling with __float128 support
CQ = gcc -DQUAD -O3 -g -Wall -Wunused -Wextra -Wno-unknown-pragmas
LMQ = -lm $(LQ)

# source code that uses GMP and MPFR
csrc            = $(filter-out \
                    ievir.c ie.c ozcrtr.c muhs.c, \
                    $(wildcard *.c))
csrc_fftw       = $(wildcard *fftw*.c)
csrc_mpfr       = $(wildcard *mpfr*.c)
csrc_gmp        = $(wildcard *gmp*.c)
csrc_gsl        = $(wildcard *gsl*.c)
csrc_extra      = $(filter-out \
		    $(csrc_fftw) $(csrc_mpfr) $(csrc_gmp) $(csrc_gsl), \
		    $(csrc))

deps = $(wildcard *.h) Makefile
bins_fftw       = $(patsubst %.c, %,   $(csrc_fftw))
bins_fftw_l     = $(patsubst %,   %_l, $(bins_fftw))
bins_fftw_q     = $(patsubst %,   %_q, $(bins_fftw))
bins_fftw_d     = $(patsubst %,   %_d, $(bins_fftw))
bins_mpfr       = $(patsubst %.c, %,   $(csrc_mpfr))
bins_mpfr_d     = $(patsubst %,   %_d, $(bins_mpfr))
bins_gmp        = $(patsubst %.c, %,   $(csrc_gmp))
bins_gmp_d      = $(patsubst %,   %_d, $(bins_gmp))
bins_gsl        = $(patsubst %.c, %,   $(csrc_gsl))
bins_gsl_l      = $(patsubst %,   %_l, $(bins_gsl))
bins_gsl_q      = $(patsubst %,   %_q, $(bins_gsl))
bins_gsl_d      = $(patsubst %,   %_d, $(bins_gsl))
bins_extra      = $(patsubst %.c, %,   $(csrc_extra))
bins_extra_d    = $(patsubst %,   %_d, $(bins_extra))
bins            = $(patsubst %.c, %,   $(csrc))
bins_d          = $(patsubst %,   %_d, $(bins))
bins_l          = $(bins_fftw_l) $(bins_gsl_l)
bins_q          = $(bins_fftw_q) $(bins_gsl_q)

# only common objects
exe: $(bins) $(bins_l) $(bins_q)

all: $(bins) $(bins_l) $(bins_q) $(bins_d)

$(bins_gsl) : %: %.c $(deps) slowdht.h
	$(CX) -DDHT=1 -DSLOWDHT=1 -o $@ $< $(LMX) $(LGSL)

$(bins_gsl_l) : %_l: %.c $(deps) slowdht.h
	$(CX) -DDHT=1 -DLDBL -DSLOWDHT=1 -o $@ $< $(LMX) $(LGSL)

$(bins_gsl_q) : %_q: %.c $(deps) slowdht.h
	$(CQ) -DDHT=1 -DF128 -DSLOWDHT=1 -o $@ $< $(LMQ) $(LGSL)

$(bins_gsl_d) : %_d: %.c $(deps) slowdht.h
	$(CC) -DDHT=1 -DSLOWDHT=1 -o $@ $< $(LM) $(LGSL)

$(bins_fftw) : %: %.c $(deps)
	$(CX) -DFFTW=1 -o $@ $< $(LMX) $(LFFTW)

$(bins_fftw_l) : %_l: %.c $(deps)
	$(CX) -DFFTW=1 -DLDBL -o $@ $< $(LMX) $(LFFTWL)

$(bins_fftw_q) : %_q: %.c $(deps)
	$(CQ) -DFFTW=1 -o $@ $< $(LQ) $(LFFTWQ)

$(bins_fftw_d) : %_d: %.c $(deps)
	$(CC) -DFFTW -o $@ $< $(LM) $(LFFTW)

$(bins_mpfr) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LMPFR)

$(bins_mpfr_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LMPFR)

$(bins_gmp) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX) $(LGMP)

$(bins_gmp_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM) $(LGMP)

$(bins_extra) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_extra_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

syncipo:
	rsync -avL *.[ch] *.py *.sh cz1@129.109.88.204:lwork/vir1

sync: syncipo
	rsync -avL *.[ch] *.py *.sh oo1@stampede.tacc.utexas.edu:/work/02464/oo1/vir

# files are kraken are automatically on keeneland
sync2::
	rsync -e gsissh -avL *.[ch] kraken-gsi.nics.xsede.org:work/vir

clean:
	rm -rf $(bins) $(bins_d) $(bins_l) $(bins_q) \
	  *~ a.out gmon.out MTSEED bin obj DHT*.dat \
	  *Bn*.dat snapshot* tmp*.dat *.pyc \
	  ieoptparam*.err ieoptparam*.out \
	  r[0-9][0-9][0-9][a-z]* vv*.dat
	rstrip.py -Rv

usb:
	$(MAKE) -C ../.. $@

pack:
	$(MAKE) -C ../.. $@

$(prj)code.zip: fft.h fftmpfr.h xdouble.h slowdht.h fftx.h zcom.h \
  		ievir.c ieutil.h \
		ievirmpfr.c ieutilmpfr.h \
		ybgvir.c kirkvir.c \
		ozcrtr.c \
		ie.c ielmv.h iemdiis.h \
		base/muhs.c \
		scifmt.py ieoptparam.py ievirextra.py \
		iebenchmark.sh \
		README.md Makefile
	zip -r $@ $^

zip: $(prj)code.zip

Bossman: zip
	rsync -avL $(prj)code.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/prog

