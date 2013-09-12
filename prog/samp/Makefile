# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm
# optimized compiling
CX = icc -g -O3 -Wall -Wunused -Wextra -march=native
LMX =

deps = dg.h dgcsep.h dgrjw.h dgdb.h zcom.h Makefile \
       mcutil.h
bins = mcrat mcgc mcgcr
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

sync::
	rsync -avL *.[ch] cz1@129.109.88.204:lwork/vir
	rsync -avL *.[ch] oo1@stampede.tacc.utexas.edu:work/vir

clean:
	rm -rf $(bins) $(bins_d) *~ a.out *.his MTSEED bin obj
