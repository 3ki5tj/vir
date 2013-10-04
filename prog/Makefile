# debugging
CC = gcc -O0 -Wall -Wunused -Wextra -g
LM = -lm
# optimized compiling
CX = icc -O3 -march=native -ipo -Wall -Wunused -Wremarks -Wcheck -wd11001 -wd11006
LMX =

deps = dg.h dgcsep.h dgsc.h dgrjw.h dgdb.h zcom.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bins)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -rf $(bins) $(bins_d) *~ a.out *.his MTSEED bin obj

usb:
	$(MAKE) -C .. usb
