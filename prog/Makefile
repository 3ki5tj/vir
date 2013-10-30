# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm
# optimized compiling
CX = icc -O3 -march=native -ipo -Wall -Wunused -Wremarks -Wcheck -wd11001 -wd11006
LMX =

deps = dg.h dgcsep.h dgsc.h dgrjw.h dgring.h dgmapl.h dgaut.h dghash.h zcom.h testutil.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

sync:
	rsync -avzL *.[ch] cz1@129.109.88.204:lwork/vir

clean:
	rm -rf $(bins) $(bins_d) *~ a.out *.his MTSEED bin obj r[0-9][0-9][0-9]hs

usb:
	$(MAKE) -C .. $@

pack:
	$(MAKE) -C .. $@

zip:
	$(MAKE) -C .. $@

