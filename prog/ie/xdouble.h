#ifndef XDOUBLE_H__
#define XDOUBLE_H__


#ifdef QUAD

#include <quadmath.h>
typedef __float128 xdouble;
#define FFTWPFX(f) fftwq_##f
#define DBLSCNF "Q"
#define DBLPRNF "Q"

#elif defined(LDBL)

typedef long double xdouble;
#define FFTWPFX(f) fftwl_##f
#define DBLSCNF "L"
#define DBLPRNF "L"

#else

typedef double xdouble;
#define FFTWPFX(f) fftw_##f
#define DBLSCNF "l"
#define DBLPRNF ""

#endif



#endif /* XDOUBLE_H__ */

