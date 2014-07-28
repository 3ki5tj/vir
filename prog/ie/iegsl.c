/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc iegsl.c -lgsl -lgslcblas
 * This program works for both even and odd dimensions
 * */
#ifndef DHT
#define DHT
#endif
#ifdef FFT0
#undef FFT0
#endif
#ifdef FFTW
#undef FFTW
#endif
#include "ievir.c"

