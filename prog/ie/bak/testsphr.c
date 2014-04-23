#include <string.h>
#include "fft.h"

/* 3D Fourier transform slow version
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0 */
static int fft3dsphr00_foo(double a[], int n, double dx, double dk, int sgn)
{
  int i, j;
  double k, x;
  double *b;

  b = malloc(n * sizeof(a[0]));
  for (i = 1; i < n; i++) {
    k = i * dk;
    b[i] = 0;
    for (j = 1; j < n; j++) {
      x = j * dx;
      b[i] += a[j] * x * sin(k * x);
    }
    b[i] *= 4 * M_PI * dx / k;
    if (sgn <= 0) b[i] /= 8 * M_PI * M_PI * M_PI;
  }
  /* skip the i == 0 term */
  for (i = 1; i < n; i++) a[i] = b[i];
  free(b);
  return 0;
}



/* 3D Fourier transform slow version
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0 */
static int fft3dsphr11_foo(double a[], int n, double dx, double dk, int sgn)
{
  int i, j;
  double k, x;
  double *b;

  b = malloc(n * sizeof(a[0]));
  for (i = 0; i < n; i++) {
    k = (i + .5) * dk;
    b[i] = 0;
    for (j = 0; j < n; j++) {
      x = (j + .5) * dx;
      b[i] += a[j] * x * sin(k * x);
    }
    b[i] *= 4 * M_PI * dx / k;
    if (sgn <= 0) b[i] /= 8 * M_PI * M_PI * M_PI;
  }
  for (i = 0; i < n; i++) a[i] = b[i];
  free(b);
  return 0;
}



#define N 8

static void cmp(double dx, double dk, const char *flag)
{
  int i;
  double a[N], a2[N];

  printf("spherical transform %s\n", flag);

  for (i = 0; i < N; i++)
    a[i] = a2[i] = (i * dx < 1) ? -1 : 0;
  if (strcmp(flag, "11") == 0) {
    fft3dsphr11(a, N, dx, dk, 1);
    fft3dsphr11_foo(a2, N, dx, dk, 1);
  } else {
    fft3dsphr00(a, N, dx, dk, 1);
    fft3dsphr00_foo(a2, N, dx, dk, 1);
  }
  for (i = 0; i < N; i++)
    printf("i %2d: %+14.6f, %14.6f\n", i, a[i], a2[i]);
  printf("\n");

  /* do the inverse transform */
  if (strcmp(flag, "11") == 0) {
    fft3dsphr11(a, N, dk, dx, -1);
    fft3dsphr11_foo(a2, N, dk, dx, -1);
  } else {
    fft3dsphr00(a, N, dk, dx, -1);
    fft3dsphr00_foo(a2, N, dk, dx, -1);
  }
  for (i = 0; i < N; i++)
    printf("i %2d: %+14.6f, %14.6f\n", i, a[i], a2[i]);
  printf("\n\n");
}



int main(void)
{
  double dx = 1/3.999, dk = M_PI / N / dx;
  cmp(dx, dk, "00");
  cmp(dx, dk, "11");
  return 0;
}


