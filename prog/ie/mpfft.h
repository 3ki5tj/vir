#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>



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
static int mpfft(mpfr_t a[], int n, int sign)
{
  mpfr_t s, c, bs, bsh, bcb, bth, tmpre, tmpim, tmp;
  int gaddr, gstep = 1; /* gaddr is the starting index of each group */
  int coupid, coups;  /* coups is how many couples in each group */
  int i, j, m;

  /* check if n is an integer power of 2 */
  for (m = n; (m & 1) == 0; m >>= 1) ;
  if (m > 1) {
    fprintf(stderr, "n %d is not a power of 2\n", n);
    exit(1);
  }

  mpfr_init(s);
  mpfr_init(c);
  mpfr_init(bs);
  mpfr_init(bsh);
  mpfr_init(bcb);
  mpfr_init(bth);
  mpfr_init(tmpre);
  mpfr_init(tmpim);
  mpfr_init(tmp);

  mpfr_const_pi(bth, MPFR_RNDN);
  if (sign > 0) mpfr_neg(bth, bth, MPFR_RNDN);

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
      mpfr_swap(a[i*2], a[j*2]);
      mpfr_swap(a[i*2+1], a[j*2+1]);
    }
  }

  /* Danielson and Lanczos */
  mpfr_set_si(bsh, 0, MPFR_RNDN);
  while (n > gstep) {
    mpfr_set(bs, bsh, MPFR_RNDN); /* bs = bsh */
    mpfr_div_si(bsh, bth, 2, MPFR_RNDN);
    mpfr_sin(bsh, bsh, MPFR_RNDN); /* bsh = sin(0.5 * bth) */
    mpfr_sqr(bcb, bsh, MPFR_RNDN);
    mpfr_mul_si(bcb, bcb, 2, MPFR_RNDN); /* bcb = 2*bsh*bsh; */
    mpfr_set_si(c, 1, MPFR_RNDN); /* c = 1; */
    mpfr_set_si(s, 0, MPFR_RNDN); /* s = 0; */
    coups = gstep;
    gstep *= 2;
    /* for each couple
     * NOTE: the value of bth (its sine and also its cosine) is only
     * relevent to the couple id, that means, two couples with a same
     * couple id will share a same theta */
    for (coupid = 0; coupid < coups; coupid++) {
      /* for each group */
      for (gaddr = 0; gaddr < n; gaddr += gstep) {
        i = gaddr + coupid; j = i + coups;
        /* calculate this couple by following fomula
         * ai = ai + aj * (c + s i);
         * aj = ai - aj * (c + s i);
         */
        mpfr_mul(tmpre, a[j*2],   c, MPFR_RNDN);
        mpfr_mul(tmp,   a[j*2+1], s, MPFR_RNDN);
        mpfr_sub(tmpre, tmpre, tmp, MPFR_RNDN); /* tmpre = re[j]*c - im[j]*s; */

        mpfr_mul(tmpim, a[j*2],   s, MPFR_RNDN);
        mpfr_mul(tmp,   a[j*2+1], c, MPFR_RNDN);
        mpfr_add(tmpim, tmpim, tmp, MPFR_RNDN); /* tmpim = re[j]*s + im[j]*c; */

        mpfr_sub(a[j*2], a[i*2], tmpre, MPFR_RNDN); /* re[j] = re[i] - tmpre; */
        mpfr_add(a[i*2], a[i*2], tmpre, MPFR_RNDN); /* re[i] += tmpre; */
        mpfr_sub(a[j*2+1], a[i*2+1], tmpim, MPFR_RNDN); /* im[j] = im[i] - tmpim; */
        mpfr_add(a[i*2+1], a[i*2+1], tmpim, MPFR_RNDN); /* im[i] += tmpim; */
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
      mpfr_set(tmp, c, MPFR_RNDN); /* tmp = c */
      mpfr_mul(tmpim, tmp, bcb, MPFR_RNDN); /* tmpim = tmp * bcb */
      mpfr_mul(tmpre, s, bs, MPFR_RNDN); /* tmpre = s * bs */
      mpfr_add(tmpim, tmpim, tmpre, MPFR_RNDN); /* tmpim = tmp * bcb + s * bc */
      mpfr_sub(c, c, tmpim, MPFR_RNDN); /* c += -(tmp = c)*bcb - s*bs; */

      mpfr_mul(tmpim, tmp, bs, MPFR_RNDN); /* tmpim = (tmp == c) * bs */
      mpfr_mul(tmpre, s, bcb, MPFR_RNDN); /* tmpre = s * bcb */
      mpfr_sub(tmpim, tmpim, tmpre, MPFR_RNDN); /* (tmp == c) * bc - s * bcb */
      mpfr_add(s, s, tmpim, MPFR_RNDN); /* s += -s*bcb + tmp*bs; */
    } /* end for each couple */
    mpfr_div_si(bth, bth, 2, MPFR_RNDN);
  }
  mpfr_clear(s);
  mpfr_clear(c);
  mpfr_clear(bs);
  mpfr_clear(bsh);
  mpfr_clear(bcb);
  mpfr_clear(bth);
  mpfr_clear(tmpre);
  mpfr_clear(tmpim);
  mpfr_clear(tmp);
  return 0;
}



mpfr_t *mpfft_arr1d_ = 0;
int mpfft_arr1d_n_ = 0;



#define MPFFT_ARR1D_FREE() \
  mpfft_arr1d_free(&mpfft_arr1d_n_, &mpfft_arr1d_)

__inline static void mpfft_arr1d_free(int *n, mpfr_t **arr)
{
  int i;

  for (i = 0; i < (*n); i++) mpfr_clear((*arr)[i]);
  if (*arr != NULL) {
    free(*arr);
    *arr = NULL;
  }
  *n = 0;
}



#define MPFFT_ARR1D_ALLOC(nnew) \
  mpfft_arr1d_alloc(nnew, &mpfft_arr1d_n_, &mpfft_arr1d_)

__inline static mpfr_t *mpfft_arr1d_alloc(int nnew,
    int *n, mpfr_t **arr)
{
  int i;

  if (nnew > (*n)) {
    for (i = 0; i < (*n); i++) mpfr_clear((*arr)[i]);
    if (*arr != NULL) free(*arr);
    *n = nnew;
    if ((*arr = malloc((*n) * sizeof(mpfr_t))) == NULL) {
      fprintf(stderr, "no memory for n %d\n", *n);
      exit(1);
    }
    for (i = 0; i < (*n); i++) mpfr_init((*arr)[i]);
  }
  return *arr;
}



/* sine transform, return \int 2 * sin(k x) a(x) dx
 * a[0] is unused and should always be zero */
__inline static int mpsint00(mpfr_t a[], int n, mpfr_t arr[])
{
  int err, i;

  if (arr == NULL) arr = MPFFT_ARR1D_ALLOC(n*4);

  for (i = 0; i < n; i++) {
    if (i != 0) {
      mpfr_set(arr[i*2], a[i], MPFR_RNDN);
      mpfr_neg(arr[(2*n - i)*2], a[i], MPFR_RNDN);
    }
    mpfr_set_si(arr[i*2+1], 0, MPFR_RNDN);
    mpfr_set_si(arr[(n + i)*2+1], 0, MPFR_RNDN);
  }
  mpfr_set_si(arr[0], 0, MPFR_RNDN);
  mpfr_set_si(arr[n*2], 0, MPFR_RNDN);

  err = mpfft(arr, n * 2, 0);

  for (i = 1; i < n; i++)
    mpfr_set(a[i], arr[i*2+1], MPFR_RNDN);
  return err;
}



/* sine transform, return \int 2 * sin(k x) a(x) dx */
__inline static int mpsint11(mpfr_t a[], int n, mpfr_t arr[])
{
  int err, i;
  mpfr_t c, s, c1, s1, th, tmp;

  mpfr_init(c);
  mpfr_init(s);
  mpfr_init(c1);
  mpfr_init(s1);
  mpfr_init(th);
  mpfr_init(tmp);
  if (arr == NULL) arr = MPFFT_ARR1D_ALLOC(n*4);

  mpfr_const_pi(th, MPFR_RNDN);
  mpfr_div_si(th, th, 2*n, MPFR_RNDN); /* th = M_PI/(2*n); */
  mpfr_cos(c1, th, MPFR_RNDN);
  mpfr_sin(s1, th, MPFR_RNDN);
  mpfr_div_si(th, th, 2, MPFR_RNDN); /* th /= 2 */
  mpfr_cos(c, th, MPFR_RNDN);
  mpfr_sin(s, th, MPFR_RNDN);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*(i+.5)/(2*n)); s = sin(M_PI*(i+.5)/(2*n)); */
    mpfr_mul(arr[i*2], a[i], c, MPFR_RNDN);
    mpfr_neg(arr[(2*n-1-i)*2], arr[i*2], MPFR_RNDN);
    mpfr_mul(arr[i*2+1], a[i], s, MPFR_RNDN);
    mpfr_set(arr[(2*n-1-i)*2+1], arr[i*2+1], MPFR_RNDN);

    /* c' = th = c * s1 - s * s1 */
    mpfr_mul(th, c, c1, MPFR_RNDN);
    mpfr_mul(tmp, s, s1, MPFR_RNDN);
    mpfr_sub(th, th, tmp, MPFR_RNDN);
    /* s = s * c1 + c * s1 */
    mpfr_mul(s, s, c1, MPFR_RNDN);
    mpfr_mul(tmp, c, s1, MPFR_RNDN);
    mpfr_add(s, s, tmp, MPFR_RNDN);
    mpfr_set(c, th, MPFR_RNDN);
  }

  err = mpfft(arr, n * 2, 0);

  mpfr_set_si(c, 1, MPFR_RNDN);
  mpfr_set_si(s, 0, MPFR_RNDN);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*i/(2*n)); s = sin(M_PI*i/(2*n)); */
    /* re' = re * c - im * s; im' = re * s + im * c; */
    mpfr_mul(th, arr[i*2], s, MPFR_RNDN);
    mpfr_mul(tmp, arr[i*2+1], c, MPFR_RNDN);
    mpfr_add(a[i], th, tmp, MPFR_RNDN);

    /* c' = th = c * s1 - s * s1 */
    mpfr_mul(th, c, c1, MPFR_RNDN);
    mpfr_mul(tmp, s, s1, MPFR_RNDN);
    mpfr_sub(th, th, tmp, MPFR_RNDN);
    /* s = s * c1 + c * s1 */
    mpfr_mul(s, s, c1, MPFR_RNDN);
    mpfr_mul(tmp, c, s1, MPFR_RNDN);
    mpfr_add(s, s, tmp, MPFR_RNDN);
    mpfr_set(c, th, MPFR_RNDN);
  }

  mpfr_clear(c);
  mpfr_clear(s);
  mpfr_clear(c1);
  mpfr_clear(s1);
  mpfr_clear(th);
  mpfr_clear(tmp);
  return err;
}



/* cosine transform, return \int 2 * cos(k x) a(x) dx */
__inline static int mpcost11(mpfr_t a[], int n, mpfr_t arr[])
{
  int err, i;
  mpfr_t c, s, c1, s1, th, tmp;

  mpfr_init(c);
  mpfr_init(s);
  mpfr_init(c1);
  mpfr_init(s1);
  mpfr_init(th);
  mpfr_init(tmp);
  if (arr == NULL) arr = MPFFT_ARR1D_ALLOC(n*4);

  mpfr_const_pi(th, MPFR_RNDN);
  mpfr_div_si(th, th, 2*n, MPFR_RNDN); /* th = M_PI/(2*n); */
  mpfr_cos(c1, th, MPFR_RNDN);
  mpfr_sin(s1, th, MPFR_RNDN);
  mpfr_div_si(th, th, 2, MPFR_RNDN); /* th /= 2 */
  mpfr_cos(c, th, MPFR_RNDN);
  mpfr_sin(s, th, MPFR_RNDN);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*(i+.5)/(2*n)); s = sin(M_PI*(i+.5)/(2*n)); */
    mpfr_mul(arr[i*2], a[i], c, MPFR_RNDN);
    mpfr_set(arr[(2*n-1-i)*2], arr[i*2], MPFR_RNDN);
    mpfr_mul(arr[i*2+1], a[i], s, MPFR_RNDN);
    mpfr_neg(arr[(2*n-1-i)*2+1], arr[i*2+1], MPFR_RNDN);

    /* c' = th = c * s1 - s * s1 */
    mpfr_mul(th, c, c1, MPFR_RNDN);
    mpfr_mul(tmp, s, s1, MPFR_RNDN);
    mpfr_sub(th, th, tmp, MPFR_RNDN);
    /* s = s * c1 + c * s1 */
    mpfr_mul(s, s, c1, MPFR_RNDN);
    mpfr_mul(tmp, c, s1, MPFR_RNDN);
    mpfr_add(s, s, tmp, MPFR_RNDN);
    mpfr_set(c, th, MPFR_RNDN);
  }

  err = mpfft(arr, n * 2, 0);

  mpfr_set_si(c, 1, MPFR_RNDN);
  mpfr_set_si(s, 0, MPFR_RNDN);
  for (i = 0; i < n; i++) {
    /* c = cos(M_PI*i/(2*n)); s = sin(M_PI*i/(2*n)); */
    /* re' = re * c - im * s; im' = re * s + im * c; */
    mpfr_mul(th, arr[i*2], c, MPFR_RNDN);
    mpfr_mul(tmp, arr[i*2+1], s, MPFR_RNDN);
    mpfr_sub(a[i], th, tmp, MPFR_RNDN);

    /* c' = th = c * s1 - s * s1 */
    mpfr_mul(th, c, c1, MPFR_RNDN);
    mpfr_mul(tmp, s, s1, MPFR_RNDN);
    mpfr_sub(th, th, tmp, MPFR_RNDN);
    /* s = s * c1 + c * s1 */
    mpfr_mul(s, s, c1, MPFR_RNDN);
    mpfr_mul(tmp, c, s1, MPFR_RNDN);
    mpfr_add(s, s, tmp, MPFR_RNDN);
    mpfr_set(c, th, MPFR_RNDN);
  }

  mpfr_clear(c);
  mpfr_clear(s);
  mpfr_clear(c1);
  mpfr_clear(s1);
  mpfr_clear(th);
  mpfr_clear(tmp);
  return err;
}



/* 3D fast Fourier transform
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0
 * a[0] is unused and should always be zero */
__inline static int mpfft3dsphr00(mpfr_t a[], int n,
    mpfr_t dx, mpfr_t dk, int sgn)
{
  int i, err = 0;
  mpfr_t fac, x;

  mpfr_init(fac);
  mpfr_init(x);

  mpfr_const_pi(x, MPFR_RNDN);
  mpfr_mul_si(x, x, 2, MPFR_RNDN); /* x = 2 * M_PI */
  if (sgn > 0) {
    mpfr_mul(fac, x, dx, MPFR_RNDN); /* fac = 2 * M_PI * dx */
  } else {
    mpfr_sqr(x, x, MPFR_RNDN);
    mpfr_div(fac, dx, x, MPFR_RNDN); /* fac = dx / (4 * M_PI * M_PI) */
  }
  mpfr_div(fac, fac, dk, MPFR_RNDN);

  for (i = 1; i < n; i++) {
    mpfr_mul_si(x, dx, i, MPFR_RNDN);
    mpfr_mul(a[i], a[i], x, MPFR_RNDN); /* form r * a(r) */
  }
  err = mpsint00(a, n, NULL); /* sine transform */
  for (i = 1; i < n; i++) {
    mpfr_div_si(x, fac, i, MPFR_RNDN);
    mpfr_mul(a[i], a[i], x, MPFR_RNDN); /* form a(k) / k */
  }
  mpfr_clear(fac);
  mpfr_clear(x);
  return err;
}



/* 3D fast Fourier transform
 * a(k) = (4 pi / k) \int {0 to +inf} r a(r) sin(k r) dr
 * for the forward transform (r --> k), set sgn to 1
 * for the inverse transform (k --> r), set sgn to -1 or 0
 * */
__inline static int mpfft3dsphr11(mpfr_t a[], int n,
    mpfr_t dx, mpfr_t dk, int sgn)
{
  int i, err = 0;
  mpfr_t fac, x;

  mpfr_init(fac);
  mpfr_init(x);

  mpfr_const_pi(x, MPFR_RNDN);
  mpfr_mul_si(x, x, 2, MPFR_RNDN); /* x = 2 * M_PI */
  if (sgn > 0) {
    mpfr_mul(fac, x, dx, MPFR_RNDN); /* fac = 2 * M_PI * dx */
  } else {
    mpfr_sqr(x, x, MPFR_RNDN);
    mpfr_div(fac, dx, x, MPFR_RNDN); /* fac = dx / (4 * M_PI * M_PI) */
  }
  mpfr_div(fac, fac, dk, MPFR_RNDN);

  for (i = 0; i < n; i++) {
    mpfr_set_si(x, i*2 + 1, MPFR_RNDN);
    mpfr_div_si(x, x, 2, MPFR_RNDN);
    mpfr_mul(x, x, dx, MPFR_RNDN); /* x = (i + 1/2) * dx */
    mpfr_mul(a[i], a[i], x, MPFR_RNDN); /* form r * a(r) */
  }
  err = mpsint11(a, n, NULL); /* sine transform */
  for (i = 0; i < n; i++) {
    mpfr_set_si(x, 2*i + 1, MPFR_RNDN);
    mpfr_div_si(x, x, 2, MPFR_RNDN);
    mpfr_div(x, fac, x, MPFR_RNDN);
    mpfr_mul(a[i], a[i], x, MPFR_RNDN); /* form a(k) / k */
  }
  mpfr_clear(fac);
  mpfr_clear(x);
  return err;
}



