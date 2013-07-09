# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm
# optimized compiling
CX = icc -O3 -Wall -Wunused -Wextra -march=native
LMX =

deps = dg.h dgsc.h dgdb.h zcom.h Makefile
bins = mcrat dgmc dgcnt
bins_d = $(patsubst %,%_d, $(bins))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bins) $(bins_d) *~ a.out *.his MTSEED

