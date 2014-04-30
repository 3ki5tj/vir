#ifndef FTED_H__
#define FTED_H__



/* even dimensional Fourier transform of
 * a spherical function using GSL */



#include <stdio.h>
#include <math.h>
#include <gsl/gsl_dht.h>



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



/* initialize dht */
__inline static gsl_dht *dht_init(int d, int Mm1, double xmax)
{
  return gsl_dht_new(Mm1, d*.5 - 1, xmax);
}


/* d-dimensional Fourier transform
 *   out(q) = \int_{0 to infinity} in(r) exp(i q r) dr
 *          = (2 pi)^(d/2) / q^(d/2-1)
 *            \int_{0 to infinity} in(r) J_{d/2-1}(q r) r^(d/2) dr
 * */
__inline static void dht_ftedsphr(gsl_dht *dht,
    double in[], double out[], int sgn, double arr[])
{
  size_t i;
  double fac, x;

  fac = pow(2*M_PI, dht->nu + 1);
  if (sgn <= 0) {
    /* gsl_dht_apply() implicitly multiplies the result by (dht->xmax)^2
     * while we need (dht->kmax)^2 */
    x = dht->kmax/dht->xmax;
    fac = (x*x)/fac;
  }
  for (i = 0; i < dht->size; i++) {
    x = (sgn > 0) ? gsl_dht_x_sample(dht, i) : gsl_dht_k_sample(dht, i);
    arr[i] = in[i] * pow(x, dht->nu);
  }
  gsl_dht_apply(dht, arr, out);
  for (i = 0; i < dht->size; i++) {
    x = (sgn > 0) ? gsl_dht_k_sample(dht, i) : gsl_dht_x_sample(dht, i);
    out[i] *= fac/pow(x, dht->nu);
  }
}



#endif

