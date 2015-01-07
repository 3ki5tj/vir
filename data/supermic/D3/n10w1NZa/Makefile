CC = icc -DDG_NORING=1 -O2 -ipo -Wall -Wremarks -wd3280 -wd11021 -openmp -g

all: a.mic a.out

deps = $(wildcard ../../*.h) Makefile

a.mic: ../../mcrat0.c $(deps)
	$(CC) -DD=3 -DN=10 -mmic $< -o $@

a.out: ../../mcrat0.c $(deps)
	$(CC) -DD=3 -DN=10       $< -o $@

