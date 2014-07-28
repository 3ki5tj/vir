/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ieodfftw.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ieodfftw.c -lfftw3l
 * Or for the 128-bit precision
 *  gcc -DF128 ieodfftw.c -lfftw3q -lquadmath -lm
 * To disable FFTW
 *  gcc -DNOFFTW ieodfftw.c -lm
 * */
#if !defined(FFTW) && !defined(FFT0)
#define FFTW
#endif
#ifdef DHT
#undef DHT
#endif
#include "ievir.c"
