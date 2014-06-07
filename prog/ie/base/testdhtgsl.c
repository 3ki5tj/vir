/* demonstration of the GSL discrete hankel transform library
 * gcc testdhtgsl.c -lgsl -lgslcblas */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
#include "fted.h"



/* brute-force implementation of the discrete Hankel transform
 *  out[m-1] = int_{0 to xmax} t dt J_{nu}(j_{nu,m} t) in(t)
 * with
 *  in[k-1] = in(j_{nu,k}/j_{nu,M})
 * the integral can be discretized as
 *  out[m-1] = 2 * (xmax/j_{nu,M})^2 sum_{k = 1 to M - 1} in[k-1]
 *             J_nu(j_{nu, m} j_{nu,k} / j_{nu,M}) / J_{nu+1}(j_{nu,k})^2
 * Cf. the third equation in Section 32.1 of the GSL manual */
static void dht_direct(double *in, double *out,
    int M, double nu, double xmax, gsl_dht *dht)
{
  int m, k;
  double x, *zeros, *jsqrs;

  zeros = calloc(M + 1, sizeof(double));
  jsqrs = calloc(M + 1, sizeof(double));

  /* compute the zeros of the Bessel function */
  for (m = 0; m <= M; m++) {
    zeros[m] = gsl_sf_bessel_zero_Jnu(nu, m);
    x = gsl_sf_bessel_Jnu(nu + 1, zeros[m]);
    jsqrs[m] = x*x;
    printf("%2d: zeros %10.6f/%10.6f jsqrs %10.6f/%10.6f\n",
        m, zeros[m], dht->j[m], jsqrs[m], dht->J2[m]);
  }

  for (m = 1; m <= M; m++) {
    out[m-1] = 0;
    for (k = 1; k <= M - 1; k++) {
      /* in[k] is evaluated at xmax * zeros[k]/zeros[M]; */
      x = gsl_sf_bessel_Jnu(nu, zeros[m]*zeros[k]/zeros[M]);
      out[m-1] += in[k-1]*x/jsqrs[k];
    }
    x = xmax/zeros[M];
    out[m-1] *= 2*x*x;
  }

  free(zeros);
  free(jsqrs);
}



/* demonstrate the discrete Hankel transform */
static void test_dht(double nu, int dht_M, double dht_xmax)
{
  gsl_dht *dht;
  double x, fac, *fr, *fk, *fr1, *fk2, *fr2;
  int i;

  printf("===== Demonstration of DFT begins =====\n");

  dht = gsl_dht_new(dht_M - 1, nu, dht_xmax);

  printf("x sampling points:\n");
  for (i = 0; i < dht_M; i++) {
    x = gsl_dht_x_sample(dht, i);
    printf("x%3d/%3d: %10.6f %10.6f\n", i, dht_M, x, x/dht_xmax);
  }
  printf("\nk sampling points:\n");
  for (i = 0; i < dht_M; i++) {
    x = gsl_dht_k_sample(dht, i);
    printf("k%3d/%3d: %10.6f %10.6f Pi\n", i, dht_M, x, x*dht_xmax/M_PI);
  }

  fr = calloc(dht_M, sizeof(double));
  fk = calloc(dht_M, sizeof(double));
  fr1 = calloc(dht_M, sizeof(double));
  fk2 = calloc(dht_M, sizeof(double));
  fr2 = calloc(dht_M, sizeof(double));

  /* prepare a function f */
  for (i = 0; i <= (int) dht->size; i++) {
    x = gsl_dht_x_sample(dht, i)/2;
    fr[i] = exp(-.5*x*x);
  }

  printf("\nThe original data\n");
  for (i = 0; i <= (int) dht->size; i++)
    printf("fr %3d/%3d: %10.6f\n", i, (int) dht->size, fr[i]);

  /* apply the discrete Hankel transform, fr --> fk */
  gsl_dht_apply(dht, fr, fk);

  /* compute dht by brute force */
  dht_direct(fr, fk2, dht_M, nu, dht_xmax, dht);

  printf("\nAfter the discrete hankel transform\n");
  for (i = 0; i <= (int) dht->size; i++)
    printf("fk %3d/%3d: %10.6f %10.6f\n",
        i, (int) dht->size, fk[i], fk2[i]);

  /* apply the inverse discrete Hankel transform */
  gsl_dht_apply(dht, fk, fr1);
  /* gsl_dht_apply() multiplies the integral by xmax^2
   * for the integral has \int_{0 to xmax} t dt ...
   * but we need \int_{0 to kmax} t dt
   * so we need to multiply a factor of kmax/xmax
   * kmax = dht->j[dht->size+1] = dht->j[dht_M] */
  fac = dht->kmax/dht_xmax;
  printf("kmax %g, xmax %g, fac %g\n", dht->kmax, dht->xmax, fac);

  /* compute the inverse dht by brute force */
  dht_direct(fk2, fr2, dht_M, nu, dht->kmax, dht);

  printf("\nAfter the inverse discrete hankel transform\n");
  for (i = 0; i <= (int) dht->size; i++) {
    fr1[i] *= fac * fac;
    printf("fr %3d/%3d: %10.6f, %10.6f vs. %10.6f\n",
        i, (int) dht->size, fr1[i], fr2[i], fr[i]);
  }

  printf("===== Demonstration of DFT ends =====\n\n\n");

  free(fr);
  free(fk);
  free(fr1);
  free(fk2);
  free(fr2);
  gsl_dht_free(dht);
}



/* brute-force implementation of the discrete Hankel transform
 *  out(q) = fac/q^{d/2-1} int_{0 to xmax} dr r^{d/2} J_{d/2-1}(q r) in(r)
 * fac = (2*pi)^(d/2) if sgn == 1, (2*pi)^(-d/2) otherwise
 * */
static void fted_direct(double *in, double *out, int sgn,
    int d, int M, double xmax)
{
  int m, k;
  double nu = d*.5-1, x, fac, *zeros, *jsqrs;

  zeros = calloc(M + 1, sizeof(double));
  jsqrs = calloc(M + 1, sizeof(double));

  /* compute the zeros of the Bessel function */
  for (m = 1; m <= M; m++) {
    zeros[m] = gsl_sf_bessel_zero_Jnu(nu, m);
    x = gsl_sf_bessel_Jnu(nu + 1, zeros[m]);
    jsqrs[m] = x*x;
    printf("%2d: zeros %10.6f jsqrs %10.6f\n", m, zeros[m], jsqrs[m]);
  }

  fac = (sgn > 0) ? pow(2*M_PI, d*.5) : pow(2*M_PI, -d*.5);
  for (m = 1; m <= M; m++) {
    out[m-1] = 0;
    for (k = 1; k <= M - 1; k++) {
      /* in[k] is evaluated at xmax * zeros[k]/zeros[M]; */
      x = gsl_sf_bessel_Jnu(nu, zeros[m]*zeros[k]/zeros[M]);
      out[m-1] += in[k-1]*pow(x, d*.5)/jsqrs[k];
    }
    x = xmax/zeros[M];
    out[m-1] *= fac*2*x*x;
  }

  free(zeros);
  free(jsqrs);
}



/* use the discrete Hankel transform to Fourier transform
 * a spherical-symmetric function */
static void test_fted(int dim, int dht_M, double dht_xmax)
{
  gsl_dht *dht;
  double x, *fr, *fk, *fr1, *fk2, *fr2, *arr;
  int i;

  printf("===== Demonstration of FTED begins =====\n");

  dht = dht_init(dim, dht_M - 1, dht_xmax);

  fr = calloc(dht_M, sizeof(double));
  fk = calloc(dht_M, sizeof(double));
  arr = calloc(dht_M, sizeof(double));
  fr1 = calloc(dht_M, sizeof(double));
  fk2 = calloc(dht_M, sizeof(double));
  fr2 = calloc(dht_M, sizeof(double));

  /* prepare a function f */
  for (i = 0; i < (int) dht->size; i++) {
    x = gsl_dht_x_sample(dht, i)/2.0;
    fr[i] = exp(-.5*x*x);
  }

  printf("\nThe original data\n");
  for (i = 0; i < (int) dht->size; i++)
    printf("fr %3d/%3d: %10.6f\n", i, (int) dht->size, fr[i]);

  /* Fourier transform, fr --> fk */
  dht_ftedsphr(dht, fr, fk, 1, arr);

  /* compute the Fourier transform by brute force */
  fted_direct(fr, fk2, 1, dim, dht_M, dht_xmax);

  printf("\nAfter the discrete hankel transform\n");
  for (i = 0; i < (int) dht->size; i++)
    printf("fk %3d/%3d: %10.6f %10.6f\n",
        i, (int) dht->size, fk[i], fk2[i]);

  /* apply the inverse Fourier transform by the discrete Hankel transform */
  dht_ftedsphr(dht, fk, fr1, -1, arr);

  /* compute the inverse Fourier transform by brute force */
  fted_direct(fk2, fr2, -1, dim, dht_M, dht->kmax);

  printf("\nAfter the inverse discrete hankel transform\n");
  for (i = 0; i < (int) dht->size; i++)
    printf("fr %3d/%3d: %10.6f, %10.6f vs. %10.6f\n",
        i, (int) dht->size, fr1[i], fr2[i], fr[i]);

  printf("===== Demonstration of FTED ends =====\n\n\n");

  free(fr);
  free(fk);
  free(arr);
  free(fr1);
  free(fk2);
  free(fr2);
  gsl_dht_free(dht);
}



int main(void)
{
  int dht_M = 8;
  double dht_xmax = 10;

  test_dht(0.5, dht_M, dht_xmax);
  test_fted(2, dht_M, dht_xmax);
  return 0;
}

