#include <string.h>
#include "mpfft.h"

/* 3D Fourier transform slow version
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0 */
static int mpfft3dsphr00_foo(mpfr_t a[], int n, mpfr_t dx, mpfr_t dk, int sgn)
{
  int i, j;
  mpfr_t k, x, kx, *b;

  mpfr_init(k);
  mpfr_init(x);
  mpfr_init(kx);
  if ((b = malloc(n * sizeof(a[0]))) == NULL) {
    fprintf(stderr, "no memory for b\n");
    exit(-1);
  }
  for (i = 0; i < n; i++) mpfr_init(b[i]);

  for (i = 1; i < n; i++) {
    mpfr_mul_si(k, dk, i, MPFR_RNDN);
    mpfr_set_si(b[i], 0, MPFR_RNDN);
    for (j = 1; j < n; j++) {
      mpfr_mul_si(x, dx, j, MPFR_RNDN); /* x = dx * j */
      mpfr_mul(kx, k, x, MPFR_RNDN); /* k * x */
      mpfr_sin(kx, kx, MPFR_RNDN); /* sin(k * x) */
      mpfr_mul(kx, kx, x, MPFR_RNDN); /* x * sin(k * x) */
      mpfr_mul(kx, kx, a[j], MPFR_RNDN);
      mpfr_add(b[i], b[i], kx, MPFR_RNDN); /* b[i] += a[j] * x * sin(k x) */
    }
    mpfr_const_pi(kx, MPFR_RNDN);
    mpfr_mul(b[i], b[i], dx, MPFR_RNDN);
    mpfr_div(b[i], b[i], k, MPFR_RNDN);
    if (sgn > 0) {
      mpfr_mul_si(b[i], b[i], 4, MPFR_RNDN);
      mpfr_mul(b[i], b[i], kx, MPFR_RNDN); /*  b[i] *= 4 * M_PI * dx / k; */
    } else {
      mpfr_sqr(kx, kx, MPFR_RNDN);
      mpfr_mul_si(kx, kx, 2, MPFR_RNDN);
      mpfr_div(b[i], b[i], kx, MPFR_RNDN); /* b[i] *= dx / k / (2 * M_PI * M_PI); */
    }
  }
  /* skip the i == 0 term */
  for (i = 1; i < n; i++)
    mpfr_set(a[i], b[i], MPFR_RNDN);
  for (i = 0; i < n; i++)
    mpfr_clear(b[i]);
  free(b);
  mpfr_clear(k);
  mpfr_clear(x);
  return 0;
}



/* 3D Fourier transform slow version
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0 */
static int mpfft3dsphr11_foo(mpfr_t a[], int n, mpfr_t dx, mpfr_t dk, int sgn)
{
  int i, j;
  mpfr_t k, x, kx, *b;

  mpfr_init(k);
  mpfr_init(x);
  mpfr_init(kx);
  if ((b = malloc(n * sizeof(a[0]))) == NULL) {
    fprintf(stderr, "no memory for b\n");
    exit(-1);
  }
  for (i = 0; i < n; i++) mpfr_init(b[i]);

  for (i = 0; i < n; i++) {
    mpfr_set_si(k, i*2 + 1, MPFR_RNDN);
    mpfr_div_si(k, k, 2, MPFR_RNDN);
    mpfr_mul(k, k, dk, MPFR_RNDN); /* k = (i + 1/2) * dk */
    mpfr_set_si(b[i], 0, MPFR_RNDN);
    for (j = 0; j < n; j++) {
      mpfr_set_si(x, j*2 + 1, MPFR_RNDN);
      mpfr_div_si(x, x, 2, MPFR_RNDN);
      mpfr_mul(x, x, dx, MPFR_RNDN); /* x = (j + 1/2) * dx */
      mpfr_mul(kx, k, x, MPFR_RNDN); /* k * x */
      mpfr_sin(kx, kx, MPFR_RNDN); /* sin(k * x) */
      mpfr_mul(kx, kx, x, MPFR_RNDN); /* x * sin(k * x) */
      mpfr_mul(kx, kx, a[j], MPFR_RNDN);
      mpfr_add(b[i], b[i], kx, MPFR_RNDN); /* b[i] += a[j] * x * sin(k x) */
    }
    mpfr_const_pi(kx, MPFR_RNDN);
    mpfr_mul(b[i], b[i], dx, MPFR_RNDN);
    mpfr_div(b[i], b[i], k, MPFR_RNDN);
    if (sgn > 0) {
      mpfr_mul_si(b[i], b[i], 4, MPFR_RNDN);
      mpfr_mul(b[i], b[i], kx, MPFR_RNDN); /*  b[i] *= 4 * M_PI * dx / k; */
    } else {
      mpfr_sqr(kx, kx, MPFR_RNDN);
      mpfr_mul_si(kx, kx, 2, MPFR_RNDN);
      mpfr_div(b[i], b[i], kx, MPFR_RNDN); /* b[i] *= dx / k / (2 * M_PI * M_PI); */
    }
  }
  for (i = 0; i < n; i++)
    mpfr_set(a[i], b[i], MPFR_RNDN);
  for (i = 0; i < n; i++)
    mpfr_clear(b[i]);
  free(b);
  mpfr_clear(k);
  mpfr_clear(x);
  return 0;
}



#define N 8

static void cmp(mpfr_t dx, mpfr_t dk, const char *flag)
{
  int i;
  mpfr_t a[N], a2[N], x;

  printf("spherical transform %s\n", flag);

  mpfr_init(x);
  for (i = 0; i < N; i++) {
    mpfr_init(a[i]);
    mpfr_init(a2[i]);
    mpfr_mul_si(x, dx, i, MPFR_RNDN);
    if (mpfr_cmp_si(x, 1) < 0) {
      mpfr_set_si(a[i], -1, MPFR_RNDN);
      mpfr_set_si(a2[i], -1, MPFR_RNDN);
    } else {
      mpfr_set_si(a[i], 0, MPFR_RNDN);
      mpfr_set_si(a2[i], 0, MPFR_RNDN);
    }
  }
  if (strcmp(flag, "00") == 0) {
    mpfft3dsphr00(a, N, dx, dk, 1);
    mpfft3dsphr00_foo(a2, N, dx, dk, 1);
  } else {
    mpfft3dsphr11(a, N, dx, dk, 1);
    mpfft3dsphr11_foo(a2, N, dx, dk, 1);
  }
  for (i = 0; i < N; i++)
    mpfr_printf("i %2d: %+14.6Rf, %14.6Rf\n", i, a[i], a2[i]);
  printf("\n");

  /* do the inverse transform */
  if (strcmp(flag, "00") == 0) {
    mpfft3dsphr00(a, N, dk, dx, -1);
    mpfft3dsphr00_foo(a2, N, dk, dx, -1);
  } else {
    mpfft3dsphr11(a, N, dk, dx, -1);
    mpfft3dsphr11_foo(a2, N, dk, dx, -1);
  }
  for (i = 0; i < N; i++)
    mpfr_printf("i %2d: %+14.6Rf, %14.6Rf\n", i, a[i], a2[i]);
  mpfr_clear(x);

  for (i = 0; i < N; i++) {
    mpfr_clear(a[i]);
    mpfr_clear(a2[i]);
  }
  printf("\n\n");
}



int main(void)
{
  mpfr_t dx, dk;

  //mpfr_set_default_prec(256);
  mpfr_init(dx);
  mpfr_set_d(dx, 1, MPFR_RNDN);
  mpfr_div_d(dx, dx, 3.999, MPFR_RNDN); /* dx = 0.25 */

  mpfr_init(dk);
  mpfr_const_pi(dk, MPFR_RNDN);
  mpfr_div(dk, dk, dx, MPFR_RNDN);
  mpfr_div_si(dk, dk, N, MPFR_RNDN); /* dk = M_PI / N / dx */

  cmp(dx, dk, "00");
  cmp(dx, dk, "11");

  mpfr_clear(dx);
  mpfr_clear(dk);
  return 0;
}


