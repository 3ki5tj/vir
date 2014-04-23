/* demonstration of the GSL discrete hankel transform library
 * gcc dht.c -lgsl -lgslcblas */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_sf_bessel_zero.h>



int dht_M = 8;
double dht_nu = 0.5;
double dht_xmax = 10;



/* brute-force implementation of the discrete Hankel transform
 * out[m-1] = 2*(xmax/j_{nu, M})^2 sum_{k = 1 to M - 1} in[k-1]
 *            J_nu(j_{nu, m} j_{nu,k} / j_{nu,M}) / J_{nu+1}(j_{nu,k})^2
 * Cf. the third equation in Section 32.1 of the GSL manual */
static void dht_foo(double *in, double *out, int M, double nu, double xmax)
{
  int m, k;
  double x, *zeros, *jsqrs;

  zeros = calloc(M + 1, sizeof(double));
  jsqrs = calloc(M + 1, sizeof(double));

  /* compute the zeros of the Bessel function */
  for (m = 1; m <= M; m++) {
    zeros[m] = gsl_sf_bessel_zero_Jnu(nu, m);
    x = gsl_sf_bessel_Jnu(nu + 1, zeros[m]);
    jsqrs[m] = x*x;
    printf("%2d: zeros %10.6f jsqrs %10.6f\n", m, zeros[m], jsqrs[m]);
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



int main(void)
{
  gsl_dht *dht;
  double x, fac, *fr, *fk, *fr2, *fk2;
  int i;

  dht = gsl_dht_new(dht_M - 1, dht_nu, dht_xmax);

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
  fr2 = calloc(dht_M, sizeof(double));
  fk2 = calloc(dht_M, sizeof(double));

  /* prepare a function f */
  for (i = 0; i < dht_M; i++) {
    x = gsl_dht_x_sample(dht, i)/2;
    fr[i] = exp(-.5*x*x);
  }

  printf("\nThe original data\n");
  for (i = 0; i < dht_M; i++) {
    printf("fr %3d/%3d: %10.6f\n", i, dht_M, fr[i]);
  }

  /* apply the discrete Hankel transform */
  gsl_dht_apply(dht, fr, fk);

  /* compute dht by brute force */
  dht_foo(fr, fk2, dht_M, dht_nu, dht_xmax);

  printf("\nAfter the discrete hankel transform\n");
  for (i = 0; i < dht_M; i++) {
    printf("fk %3d/%3d: %10.6f %10.6f\n", i, dht_M, fk[i], fk2[i]);
  }

  /* apply the inverse discrete Hankel transform */
  gsl_dht_apply(dht, fk, fr2);
  fac = gsl_dht_k_sample(dht, dht_M - 1)/dht_xmax;

  printf("\nAfter the inverse discrete hankel transform\n");
  for (i = 0; i < dht_M; i++) {
    fr2[i] *= fac * fac;
    printf("fr %3d/%3d: %10.6f vs. %10.6f\n", i, dht_M, fr2[i], fr[i]);
  }

  free(fr);
  free(fk);
  free(fr2);
  free(fk2);
  gsl_dht_free(dht);
  return 0;
}
