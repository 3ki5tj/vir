#ifndef SLOWDHT_H__
#define SLOWDHT_H__



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xdouble.h"
#include <gsl/gsl_sf_bessel.h>



#ifdef SLOWDHT

#define xdht slowdht
#define XDHT(f) slowdht_ ## f

#else

#include <gsl/gsl_dht.h>
#define xdht gsl_dht
#define XDHT(f) gsl_dht_ ## f

#endif



typedef struct slowdht_struct {
  size_t    size;
  xdouble   nu;
  xdouble   inu;
  xdouble   xmax;
  xdouble   kmax;
  xdouble   *j;    /* j_{nu, s} = j[s] */
  xdouble   *J2;   /* J_{nu+1}^2(j_m) */
} slowdht;



/* J(nu, x) */
__inline xdouble besselJnu(xdouble nu, xdouble x)
{
  int inu = (int)(nu + .5);
  if (FABS(nu - inu) < 1e-7)
    return JN(inu, x);
  else
    return gsl_sf_bessel_Jnu(nu, x);
}



/* return dJ(nu, x) / dx */
__inline xdouble besseldJnu(xdouble nu, xdouble x)
{
  int inu = (int) (nu + .5);
  xdouble x1, x2;
  if (FABS(nu - inu) < 1e-7) {
    if (inu == 0) return -J1(x);
    x1 = JN(inu - 1, x);
    x2 = JN(inu + 1, x);
  } else {
    x1 = gsl_sf_bessel_Jnu(nu - 1, x);
    x2 = gsl_sf_bessel_Jnu(nu + 1, x);
  }
  return (x1 - x2) / 2;
}



/* compute the mth zero of the Bessel function
 * for xdouble == long double or __float128
 * use the GSL function get an estimate,
 * then refine it by Newton-Raphson's method  */
__inline xdouble besselzeroJnu(xdouble nu, int m)
{
  xdouble x, x1, y, yabs, dy, y1, y1abs;
  int i;

  x = gsl_sf_bessel_zero_Jnu(nu, m);
  /* the Newton-Raphson's method fails at x = 0 */
  if (FABS(x) < 1e-8) return x;
  y = besselJnu(nu, x);
  yabs = FABS(y);
  //printf("x %.20" XDBLPRNF "f, y %g, nu %g\n", x, (double) y, (double) nu);
  for ( i = 0; ; i++ ) {
    dy = besseldJnu(nu, x);
    x1 = x - y/dy;
    y1 = besselJnu(nu, x1);
    y1abs = FABS(y1);
    if (y1abs >= yabs) break;
    x = x1;
    y = y1;
    yabs = y1abs;
  }
  //printf("iter %d, x %.20" XDBLPRNF "f, y %g\n", i, x, (double) y); getchar();
  return x;
}



__inline slowdht *slowdht_new(size_t size, xdouble nu, xdouble xmax)
{
  slowdht *dht;
  size_t i;
  xdouble x;

  if ((dht = calloc(1, sizeof(*dht))) == NULL) {
    fprintf(stderr, "no memory for dht\n");
    exit(1);
  }
  dht->size = size;
  dht->nu = nu;
  dht->inu = (int) (nu + .5);
  if (FABS(dht->nu - dht->inu) > 1e-3)
    dht->inu = -1;
  dht->xmax = xmax;

  if ((dht->j = calloc(dht->size + 2, sizeof(xdouble))) == NULL) {
    fprintf(stderr, "no memory for dht->j %d\n", (int) dht->size);
    exit(1);
  }

  /* compute the zeros of the Bessel functions
   * Note, the 0th zero of nu = 0 does not exist */
  dht->j[0] = 0;
  for ( i = 1; i <= size + 1; i++ ) {
    dht->j[i] = besselzeroJnu(nu, i);
  }
  dht->kmax = dht->j[dht->size+1] / xmax;

  if ((dht->J2 = calloc(dht->size + 1, sizeof(xdouble))) == NULL) {
    fprintf(stderr, "no memory for dht->J2 %d\n", (int) dht->size);
    exit(1);
  }
  for ( i = 0; i <= size; i++ ) {
    x = besselJnu(nu + 1, dht->j[i]);
    dht->J2[i] = x*x;
  }
  return dht;
}



__inline xdouble slowdht_x_sample(slowdht *dht, int n)
{
  if ( (size_t) n <= dht->size )
    return dht->xmax * dht->j[n+1] / dht->j[dht->size+1];
  return 0;
}



__inline xdouble slowdht_k_sample(slowdht *dht, int n)
{
  if ( (size_t) n <= dht->size )
    return dht->j[n+1] / dht->xmax;
  return 0;
}



__inline int slowdht_apply(const slowdht *dht, xdouble *inp, xdouble *out)
{
  int m, k, M = (int) dht->size + 1;
  xdouble x, y;

  for ( m = 1; m < M; m++ ) {
    y = 0;
    for ( k = 1;  k < M; k++ ) {
      x = dht->j[m] * dht->j[k] / dht->j[M];
      x = besselJnu(dht->nu, x);
      y += inp[k - 1] * x / dht->J2[k];
    }
    x = dht->xmax / dht->j[M];
    out[m - 1] = 2*(x*x)*y;
  }
  return 0;
}



__inline void slowdht_free(slowdht *dht)
{
  free(dht->j);
  free(dht->J2);
  free(dht);
}



#endif /* SLOW_DHT_H__ */

