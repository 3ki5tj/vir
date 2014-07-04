#ifndef XDOUBLE_H__
#define XDOUBLE_H__


#if (defined(QUAD) || defined(F128))

/* support __float128 for GCC
 * Intel compiler defines __GNUC__, but it does not have __float128 */
#if defined(__GNUC__) && !defined(__INTEL_COMPILER) \
  && ( (__GNUC__ == 4 && __GNUC_MINOR__ >= 6) || __GNUC__ > 4 )
  #ifndef HAVEF128
  #define HAVEF128 1
  #endif
  #include <quadmath.h>
  /* ignore warnings for "%Qf"
   * assume `diagnostic push' is available in this case */
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wformat"
  #pragma GCC diagnostic ignored "-Wformat-extra-args"
#else
  #ifndef HAVEF128
  #define HAVEF128 0
  #endif
  #error "Your compiler does not support __float128, if it does, please change xdouble.h."
#endif


#include <quadmath.h>
typedef __float128 xdouble;
#define FFTWPFX(f) fftwq_##f
#define XDBLSCNF "Q"
#define XDBLPRNF "Q"
#define XDBLPRNE "%41.32Qe" /* full precision print */
#define STRPREC "f128"
#define PI M_PIq
#define SQRT(x)   sqrtq(x)
#define EXP(x)    expq(x)
#define SIN(x)    sinq(x)
#define COS(x)    cosq(x)
#define LOG(x)    logq(x)
#define POW(x, y) powq(x, y)
#define J0(x)     j0q(x)
#define J1(x)     j1q(x)
#define JN(n, x)  jnq(n, x)



#elif (defined(LDBL) || defined(LONG))

#ifndef HAVEF128
#define HAVEF128 0
#endif

typedef long double xdouble;
#define FFTWPFX(f) fftwl_##f
#define XDBLSCNF "L"
#define XDBLPRNF "L"
#define XDBLPRNE "%27.18Le"
#define STRPREC "ldbl"
#define PI (xdouble) 3.1415926535897932384626433832795L
#define SQRT(x)   sqrtl(x)
#define EXP(x)    expl(x)
#define SIN(x)    sinl(x)
#define COS(x)    cosl(x)
#define LOG(x)    logl(x)
#define POW(x, y) powl(x, y)
#define J0(x)     j0l(x)
#define J1(x)     j1l(x)
#define JN(n, x)  jnl(n, x)

#else

#ifndef HAVEF128
#define HAVEF128 0
#endif

typedef double xdouble;
#define FFTWPFX(f) fftw_##f
#define XDBLSCNF "l"
#define XDBLPRNF ""
#define XDBLPRNE "%22.14e"
#define STRPREC ""
#define PI 3.1415926535897932384626433832795
#define SQRT(x)   sqrt(x)
#define EXP(x)    exp(x)
#define SIN(x)    sin(x)
#define COS(x)    cos(x)
#define LOG(x)    log(x)
#define POW(x, y) pow(x, y)
#define J0(x)     j0(x)
#define J1(x)     j1(x)
#define JN(n, x)  jn(n, x)

#endif

__inline static xdouble FABS(xdouble x)
{
  return x < 0 ? -x : x;
}

__inline static xdouble pow_si(xdouble x, int n)
{
  int sgn = 1;
  xdouble y;

  if ( n < 0) { sgn = -1; n = -n; }
  for ( y = 1; n; n-- ) y *= x;
  return sgn > 0 ? y : 1/y;
}



#if HAVEF128
  #pragma GCC diagnostic pop
#endif



#endif /* XDOUBLE_H__ */

