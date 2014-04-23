#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



/* Discription:
 *   replace the complex array a[] by their discrete fast fourier transform
 *
 * Parameter:
 *   o a[] is both the input and the output
 *   o a[i*2] and a[i*2+1] are the real and complex parts, respectively
 *   o n must be a power of 2
 *   o if sign is 1, do fast fourier transform,
 *     if sign is 0 or -1, do the inverse transform
 *
 * Return:
 *   on sucess, returns 0;
 *   if n isn't an integer power of 2, returns 1
 */
static int fft(double a[], int n, int sign)
{
  double s, c, bs, bsh, bcb, bth = M_PI, tmpre, tmpim, tmp;
  int gaddr, gspan = 1; /* gaddr is the starting index of each group */
  int coupid, coups;  /* coups is how many couples in each group */
  int i, j, m;

  /* check if n is an integer power of 2 */
  for (m = n; (m & 1) == 0; m >>= 1) ;
  if (m > 1) {
    fprintf(stderr, "n %d is not a power of 2\n", n);
    exit(1);
  }

  if (sign > 0) bth = -bth;

  j = 0;
  /* bit reversal */
  for (i = 1; i < n - 1; i++) {
    m = n >> 1;
    while (j >= m) {
      j -= m;
      m >>= 1;
    }
    j += m;
    if (j > i) {
      tmp = a[i*2]; a[i*2] = a[j*2]; a[j*2] = tmp;
      tmp = a[i*2+1]; a[i*2+1] = a[j*2+1]; a[j*2+1] = tmp;
    }
  }

  /* Danielson and Lanczos */
  bsh = 0;
  while (n > gspan)
  {
    bs = bsh;
    bsh = sin(0.5 * bth);
    bcb = 2*bsh*bsh;
    c = 1;
    s = 0;
    coups = gspan;
    gspan *= 2;
    /* for each couple
     * NOTE: the value of bth (its sine and also its cosine) is only
     * relevent to the couple id, that means, two couples with a same
     * couple id will share a same theta */
    for (coupid = 0; coupid < coups; coupid++) {
      /* for each group */
      for (gaddr = 0; gaddr < n; gaddr += gspan) {
        i = gaddr + coupid; j = i + coups;
        /* calculate this couple by following fomula
         * ai = ai + aj * (c + s i);
         * aj = ai - aj * (c + s i);
         */
        tmpre = a[j*2]*c - a[j*2+1]*s;
        tmpim = a[j*2]*s + a[j*2+1]*c;
        a[j*2] = a[i*2] - tmpre;
        a[i*2] += tmpre;
        a[j*2+1] = a[i*2+1] - tmpim;
        a[i*2+1] += tmpim;
      }
      /* deal next couple, increase  cosine and sine by
       * cos(a+b) = cos(a)*cos(b) - sin(a)*sin(b);
       * sin(a+b) = sin(a)*cos(b) + cos(a)*sin(b);
       * but when b is very small(when n is a verylarge number),
       * cos(b) will be very close to 1 and inaccute,
       * so we replace these fomulas by introduce bcb = 1-cos(b) = 2*sin(b/2)*sin(b/2),
       * -- a reletively small number,
       * cos(a+b) = cos(a) + (cos(a)*(-bcb) - sin(a)*sin(b));
       * sin(a+b) = sin(a) + (sin(a)*(-bcb) + cos(a)*sin(b));
       */
      c += -(tmp = c)*bcb - s*bs;
      s += -s*bcb + tmp*bs;
    } /* end for each couple */
    bth /= 2;
  }
  return 0;
}



/* sine transform, return \int 2 * sin(k x) a(x) dx
 * a[0] is unused and should always be zero */
__inline static int sint00(double a[], int n)
{
  int err, i;
  double *arr;

  if ((arr = malloc(n * 4 * sizeof(*a))) == NULL) {
    fprintf(stderr, "no memory for n %d\n", n);
    exit(1);
  }
  for (i = 0; i < n; i++) {
    if (i != 0) {
      arr[i*2] = a[i];
      arr[(2*n - i)*2] = -a[i];
    }
    arr[(n + i)*2+1] = arr[i*2+1] = 0;
  }
  arr[0] = arr[n*2] = 0;
  err = fft(arr, n*2, 0);
  for (i = 1; i < n; i++) a[i] = arr[i*2+1];
  free(arr);
  return err;
}



/* sine transform, return \int 2 * sin(k x) a(x) dx */
__inline static int sint11(double a[], int n)
{
  int err, i;
  double *arr, c, s, th, c1, s1, tmp;

  if ((arr = malloc(n * 4 * sizeof(*a))) == NULL) {
    fprintf(stderr, "no memory for n %d\n", n);
    exit(1);
  }

  th = M_PI/(2*n);
  c1 = cos(th);
  s1 = sin(th);
  th /= 2;
  c = cos(th);
  s = sin(th);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*(i+.5)/(2*n)); s = sin(M_PI*(i+.5)/(2*n)); */
    arr[i*2] = a[i] * c;
    arr[(2*n-1-i)*2] = -arr[i*2];
    arr[i*2+1] = a[i] * s;
    arr[(2*n-1-i)*2+1] = arr[i*2+1];
    tmp = c * c1 - s * s1;
    s = s * c1 + c * s1;
    c = tmp;
  }
  err = fft(arr, 2 * n, 0);
  c = 1;
  s = 0;
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*i/(2*n)); s = sin(M_PI*i/(2*n)); */
    /* re' = re * c - im * s; im' = re * s + im * c; */
    a[i] = arr[i*2] * s + arr[i*2+1] * c;
    tmp = c * c1 - s * s1;
    s = s * c1 + c * s1;
    c = tmp;
  }
  free(arr);
  return err;
}



/* cosine transform, return \int 2 * cos(k x) a(x) dx */
__inline static int cost11(double a[], int n)
{
  int err, i;
  double *arr, c, s, th, c1, s1, tmp;

  if ((arr = malloc(n * 4 * sizeof(*a))) == NULL) {
    fprintf(stderr, "no memory for n %d\n", n);
    exit(1);
  }

  th = M_PI/(2*n);
  c1 = cos(th);
  s1 = sin(th);
  th /= 2;
  c = cos(th);
  s = sin(th);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*(i+.5)/(2*n)); s = sin(M_PI*(i+.5)/(2*n)); */
    arr[i*2] = a[i] * c;
    arr[(2*n-1-i)*2] = arr[i*2];
    arr[i*2+1] = a[i] * s;
    arr[(2*n-1-i)*2+1] = -arr[i*2+1];
    tmp = c * c1 - s * s1;
    s = s * c1 + c * s1;
    c = tmp;
  }
  err = fft(arr, 2 * n, 0);
  c = 1;
  s = 0;
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*i/(2*n)); s = sin(M_PI*i/(2*n)); */
    /* re' = re * c - im * s; im' = re * s + im * c; */
    a[i] = arr[i*2] * c - arr[i*2+1] * s;
    tmp = c * c1 - s * s1;
    s = s * c1 + c * s1;
    c = tmp;
  }
  free(arr);
  return err;
}



/* 3D fast Fourier transform
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0
 * a[0] is unused and should always be zero */
__inline static int fft3dsphr00(double a[], int n,
    double dx, double dk, int sgn)
{
  int i, err;
  double fac = 2 * M_PI * dx, x;

  if (sgn <= 0) fac /= 8 * M_PI * M_PI * M_PI;
  for (i = 1; i < n; i++)
    a[i] *= x = dx * i; /* form r * a(r) */
  err = sint00(a, n);
  for (i = 1; i < n; i++)
    a[i] *= fac / (dk * i); /* form a(k) / k */
  return err;
}



/* 3D fast Fourier transform
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0 */
__inline static int fft3dsphr11(double a[], int n,
    double dx, double dk, int sgn)
{
  int i, err;
  double fac = 2 * M_PI * dx, x;

  if (sgn <= 0) fac /= 8 * M_PI * M_PI * M_PI;
  for (i = 0; i < n; i++)
    a[i] *= x = dx * (i + .5); /* form r * a(r) */
  err = sint11(a, n);
  for (i = 0; i < n; i++)
    a[i] *= fac / (dk * (i + .5)); /* form a(k) / k */
  return err;
}
