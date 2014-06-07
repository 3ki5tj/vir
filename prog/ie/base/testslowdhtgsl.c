/* comparison of gsl_dht and slowdht
 * gcc testslowdhtgsl.c -lgsl -lgslcbas */
#include "slowdht.h"
#include <gsl/gsl_dht.h>



static void test_slowdht(double nu, int M, double xmax)
{
  gsl_dht *t;
  slowdht *s;
  double xt, *frt, *fkt;
  xdouble xs, *frs, *fks;
  int i;

  t = gsl_dht_new(M - 1, nu, xmax);
  s = slowdht_new(M - 1, nu, xmax);

  frt = calloc(M, sizeof(*frt));
  frs = calloc(M, sizeof(*frs));
  fkt = calloc(M, sizeof(*frt));
  fks = calloc(M, sizeof(*frs));

  /* compare sampling points */
  printf("x sampling points:\n");
  for (i = 0; i < M; i++) {
    xt = gsl_dht_x_sample(t, i);
    xs = slowdht_x_sample(s, i);
    printf("x%3d/%3d: %10.6f %10.6f | %10.6f %10.6f\n",
        i, M, xt, (double) xs, t->j[i], (double) s->j[i]);
  }

  printf("\nk sampling points:\n");
  for (i = 0; i < M; i++) {
    xt = gsl_dht_k_sample(t, i);
    xs = slowdht_k_sample(s, i);
    printf("k%3d/%3d: %10.6f %10.6f\n",
        i, M, xt, (double) xs);
  }

  printf("\nJ2 comparison:\n");
  for (i = 0; i < M; i++) {
    printf("J2%3d/%3d: %10.6f %10.6f\n",
        i, M, t->J2[i], (double) t->J2[i]);
  }

  /* prepare a function f */
  for (i = 0; i < M; i++) {
    xt = gsl_dht_x_sample(t, i)/2;
    frt[i] = exp(-xt*xt/2);
    xs = slowdht_x_sample(s, i)/2;
    frs[i] = EXP(-xs*xs/2);
  }

  printf("\nThe original data\n");
  for (i = 0; i < M; i++)
    printf("fr %3d/%3d: %10.6f %10.6f\n",
        i, M, frt[i], (double) frs[i]);

  /* apply the discrete Hankel transform, fr --> fk */
  gsl_dht_apply(t, frt, fkt);
  slowdht_apply(s, frs, fks);

  printf("\nAfter the discrete Hankel transform\n");
  for (i = 0; i < M; i++)
    printf("fk %3d/%3d: %10.6f %10.6f\n",
        i, M, fkt[i], (double) fks[i]);

  gsl_dht_free(t);
  slowdht_free(s);
  free(frt);
  free(frs);
  free(fkt);
  free(fks);
}



int main(void)
{
  int dht_M = 8;
  double dht_xmax = 10;

  test_slowdht(1, dht_M, dht_xmax);
  return 0;
}

