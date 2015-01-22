gcc -O0 -g -Wall -DLDBL ybgvir.c -lfftw3l -lm -o ybgvir_l
gcc -O0 -g -Wall -DLDBL kirkvir.c -lfftw3l -lm -o kirkvir_l
gcc -O0 -g -Wall -DLDBL -DFFTW ievir.c -lfftw3l -lm -o ievirfftw_l
