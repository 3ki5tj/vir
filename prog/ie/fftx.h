#ifndef FFTX_H__
#define FFTX_H__



/* convenient wrapper for high-dimensioanl Fourier transform
 * of a spherical-symmetric function
 *
 * This header uses
 *   fft.h
 *   slowdht.h
 *   ieutil.h
 *   xdouble.h (included in the above)
 * */



#include "xdouble.h"
#include <stdlib.h>



#if !defined(FFT0) && !defined(FFTW) && !defined(DHT)
#define FFT0 /* default branch, no external library */
#endif

#if defined(FFT0) || defined(FFTW)
#define FFT
#endif

#ifdef FFT0  /* home-made FFT */
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#endif

#ifdef FFTW /* for odd dimensions */
#include <fftw3.h>
#endif

#ifdef DHT /* for even and odd dimensions */
#include "slowdht.h"
int dhtdisk = SLOWDHT_USEDISK;
typedef void *FFTWPFX(plan);
#endif



#ifndef xnew
#define xnew(x, n) { \
  if ((x = calloc(n, sizeof(*(x)))) == NULL) { \
    fprintf(stderr, "no memory for %s x %d\n", #x, (int) (n)); \
    exit(1); } }
#endif

#ifndef MAKE1DARR
#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (int) (n); i_++) arr[i_] = 0; }
#endif

#ifndef FREE1DARR
#define FREE1DARR(arr, n) if ( (arr) != NULL ) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) (arr)[i_] = 0; \
  free(arr); }
#endif

#ifndef MAKE2DARR
#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }
#endif

#ifndef FREE2DARR
#define FREE2DARR(arr, n1, n2) if ( (arr) != NULL ) { \
  FREE1DARR((arr)[0], (n1) * (n2)); \
  free(arr); }
#endif



typedef struct {
  int dim;
  int K; /* K = (dim - 1)/2 for an odd dimension */
  int npt; /* number of points along r */
  int ffttype; /* 00 or 11, if sampling half grid places */
  xdouble surfr, surfk;
  xdouble facr2k, fack2r;
  xdouble B2hs, B2g;
  xdouble rmax, dr, dk;
  int dm; /* boundary of the hard sphere */
  xdouble *arr;
  xdouble *ri;
  xdouble *ki;
  xdouble *rDm1;
  xdouble *kDm1;

  /* the next four arrays are only for FFT */
  xdouble **r2pow;
  xdouble **k2pow;
  xdouble **invr2pow;
  xdouble **invk2pow;
  long *coef; /* coefficients of the spherical bessel functions */
  FFTWPFX(plan) plans[2];

  xdouble *r2p, *k2p;
#ifdef DHT
  xdht *dht;
#endif
} sphr_t;



/*
 * Interface:
 *    sphr = sphr_open(dim, npt, rmax, Rmax, ffttype);
 *    sphr_r2k(sphr, fr, fk);
 *    sphr_k2r(sphr, fk, fr);
 *    sphr_close(sphr);
 **/



/* low-level initialization
 * use sphr_open() instead */
__inline static sphr_t *sphr_init(int dim, int npt)
{
  int i;
  sphr_t *sphr;

  xnew(sphr, 1);

  sphr->dim = dim;
  sphr->K = (dim - 1)/2;
  sphr->npt = npt;
  sphr->dm = 0;

  /* B2 = (1/2) (PI*2)^(dim/2) / dim!! for an even dim
   *    = (PI*2)^((dim - 1)/2) / dim!! for an odd dim */
  sphr->B2hs = dim % 2 ? 1 : 0.5;
  for ( i = 2 + dim % 2; i <= dim; i += 2 )
    sphr->B2hs *= PI*2/i;

  /* B2 of the Gaussian model */
  sphr->B2g = SQRT(pow_si(PI, dim))/2;

  sphr->surfr = sphr->B2hs * 2 * dim;
  sphr->surfk = sphr->surfr / pow_si(PI*2, dim);

  /* auxiliary array for FFTW or DHT
   * needs npt + 1 elements for FFTW_REDFT00 */
  MAKE1DARR(sphr->arr,  npt + 1);

  MAKE1DARR(sphr->ri,   npt);
  MAKE1DARR(sphr->ki,   npt);
  MAKE1DARR(sphr->rDm1, npt);
  MAKE1DARR(sphr->kDm1, npt);

  return sphr;
}



__inline static void sphr_close(sphr_t *sphr)
{
  int npt = sphr->npt, K = sphr->K;

  FREE1DARR(sphr->arr,  npt + 1);
  FREE1DARR(sphr->ri,   npt);
  FREE1DARR(sphr->ki,   npt);
  FREE1DARR(sphr->rDm1, npt);
  FREE1DARR(sphr->kDm1, npt);
  FREE2DARR(sphr->r2pow, K, npt);
  FREE2DARR(sphr->k2pow, K, npt);
  FREE2DARR(sphr->invr2pow, K, npt);
  FREE2DARR(sphr->invk2pow, K, npt);
  free(sphr->coef);

  FREE1DARR(sphr->r2p,  npt);
  FREE1DARR(sphr->k2p,  npt);
#ifdef DHT
  XDHT(free)(sphr->dht);
#endif

#ifdef FFTW
  FFTWPFX(destroy_plan)(sphr->plans[0]);
  FFTWPFX(destroy_plan)(sphr->plans[1]);
  FFTWPFX(cleanup)();
#endif
  free(sphr);
}



/* for both FFT0 and FFTW, only for odd dimensions */
#ifdef FFT



/* get the coefficient c_l of spherical Bessel function jn(x)
 *  jn(x) = Sum_{l = 0 to n} c_l [sin(x) or cos(x)] / x^{n + l + 1},
 * where the l'th term is sin(x) if n + l is even, or cos(x) otherwise */
static void sphr_getjn(long *c, int n)
{
  int i, k;
  const char *fs[2] = {"sin(x)", "cos(x)"};

  c[0] = 1; /* j0 = sin(x)/x; */
  if ( n == -1 ) c[0] = 2; /* if n == -1 (dim == 1), 2 cos(x) */
  /* j_n(x)/x^n = (-1/x d/dx)^n j_0(x) */
  for ( k = 1; k <= n; k++ ) { /* k'th round */
    c[k] = 0;
    for ( i = k - 1; i >= 0; i-- ) {
      c[i + 1] += c[i] * (k + i);
      c[i] *= 1 - (k + i) % 2 * 2;
    }
  }
  printf("j%d(x) =", n);
  for ( i = 0; i <= n; i++ )
    printf(" %+ld*%s/x^%d", c[i], fs[(i+n+2)%2], n + i + 1);
  printf("\n");
}



/* compute
 *    out(k) = 2 fac Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' are used in the intermediate steps by the FFTW plans `p'
 * */
static void sphr_fft(int K, int npt, xdouble *in, xdouble *out,
    xdouble fac, FFTWPFX(plan) p[2], xdouble *arr, long *coef,
    xdouble **r2p, xdouble **k2q, int ffttype)
{
  int i, l, iscos;

  (void) p;

  if ( ffttype != 0 && ffttype != 1 ) return;

  /* clear the output */
  for ( i = 0; i < npt; i++ ) out[i] = 0;

  /* several rounds of transforms */
  for ( l = 0; l < K || l == 0; l++ ) {
    /* decide if we want to do a sine or cosine transform */
    iscos = (K + l + 1) % 2;

    for ( i = 0; i < npt; i++ ) {
      /* arr[i] = in[i] * pow(dx*(2*i + 1)/2, K - l) * coef[l]; */
      arr[i] = in[i] * ( K > 0 ? coef[l] * r2p[l][i] : 1);
    }

    /* NOTE: a factor of two is included in the cosine/sine transform */
#ifdef FFTW /* use FFTW */
    FFTWPFX(execute)(p[iscos]);
#else /* use the home-made FFT */
    if ( iscos ) {
      if ( ffttype == 1 ) {
        cost11(arr, npt);
      } else {
        cost00(arr, npt);
      }
    } else {
      if ( ffttype == 1 ) {
        sint11(arr, npt);
      } else {
        sint00(arr, npt);
      }
    }
#endif /* !defined(FFTW) */

    if ( iscos && ffttype == 0 ) arr[npt] = 0;
    for ( i = 0; i < npt; i++ ) {
      /* out[i] += arr[i] * fac / pow(dk*(2*i + 1)/2, K + l); */
      out[i] += arr[i] * fac * (K > 0 ? k2q[l][i] : 1);
    }
  }
}



/* convenience macros */
#define sphr_r2k(sphr, in, out) sphr_fft(sphr->K, sphr->npt, \
    in, out, sphr->facr2k, sphr->plans, sphr->arr, \
    sphr->coef, sphr->r2pow, sphr->invk2pow, sphr->ffttype)

#define sphr_k2r(sphr, in, out) sphr_fft(sphr->K, sphr->npt, \
    in, out, sphr->fack2r, sphr->plans, sphr->arr, \
    sphr->coef, sphr->k2pow, sphr->invr2pow, sphr->ffttype)



/* adjust sphr->rmax such that r = 1 lies in precision
 * at the middle point of two bins */
__inline static void sphr_adjustrmax(sphr_t *sphr, xdouble rmax)
{
  /* fix dr such that r = 1 at the middle of bins dm - 1 and dm */
  if ( sphr->ffttype ) { /* r(i) = dr * (i + .5) */
    sphr->dm = (int) (sphr->npt / rmax + 1e-8); /* number of bins in the core */
    sphr->dr = (xdouble) 1 / sphr->dm;
  } else {
    sphr->dm = (int) (sphr->npt / rmax + .5 + 1e-8);
    sphr->dr = (xdouble) 2 / (sphr->dm*2 - 1);
  }
  sphr->rmax = sphr->dr * sphr->npt;
  printf("%d bins, %d within the hard core (dr %g), rmax %g\n",
      sphr->npt, sphr->dm, (double) sphr->dr, (double) sphr->rmax);
}



#define sphr_open(dim, npt, rmax, Rmax, ffttype) \
  sphr_openfft(dim, npt, rmax, ffttype)

/* prepare data for Fourier transform of a spherical function */
__inline static sphr_t *sphr_openfft(int dim, int npt, xdouble rmax, int ffttype)
{
  int i, l, K = (dim - 1)/2;
  xdouble rl, invrl, kl, invkl;
  sphr_t *sphr;

  if ( dim % 2 != 1 ) {
    fprintf(stderr, "cannot use FFT for D = %d\n", dim);
    exit(1);
  }

  sphr = sphr_init(dim, npt);
  sphr->ffttype = ffttype;

  sphr_adjustrmax(sphr, rmax);
  sphr->dk = PI / sphr->dr / npt;

  /* do this after dr and dk are known */
  sphr->facr2k = pow_si(PI*2, K) * sphr->dr;
  sphr->fack2r = pow_si(PI*2, -K-1) * sphr->dk;

  if ( K > 0 ) { /* to avoid the D = 1 case */
    MAKE2DARR(sphr->r2pow, K, npt) /* r^{K - l} for l = 0, ..., K */
    MAKE2DARR(sphr->k2pow, K, npt) /* k^{K - l} for l = 0, ..., K */
    MAKE2DARR(sphr->invr2pow, K, npt) /* r^{-K-l} for l = 0, ..., K */
    MAKE2DARR(sphr->invk2pow, K, npt) /* k^{-K-l} for l = 0, ..., K */
  }
  for ( i = 0; i < npt; i++ ) {
    sphr->ri[i]  = sphr->dr * (i*2 + (ffttype ? 1 : 0))/2;
    sphr->ki[i]  = sphr->dk * (i*2 + (ffttype ? 1 : 0))/2;
    rl = 1;
    kl = 1;
    for ( l = 1; l <= K; l++ ) {
      sphr->r2pow[K - l][i] = (rl *= sphr->ri[i]);
      sphr->k2pow[K - l][i] = (kl *= sphr->ki[i]);
    }

    if ( ffttype == 0 && i == 0 ) continue;

    invrl = 1/rl;
    invkl = 1/kl;
    for ( l = 0; l < K; l++ ) {
      sphr->invr2pow[l][i] = invrl;
      sphr->invk2pow[l][i] = invkl;
      invrl /= sphr->ri[i];
      invkl /= sphr->ki[i];
    }

    sphr->rDm1[i] = sphr->surfr * pow_si(sphr->ri[i], dim - 1) * sphr->dr;
    sphr->kDm1[i] = sphr->surfk * pow_si(sphr->ki[i], dim - 1) * sphr->dk;
  }

  if ( K > 0 ) {
    xnew(sphr->coef, K);
    sphr_getjn(sphr->coef, K - 1);
  }

#ifdef FFTW
  /* plans[0] is the sine transform, plans[1] is the cosine transform */
  if ( ffttype ) {
    sphr->plans[0] = FFTWPFX(plan_r2r_1d)(npt, sphr->arr, sphr->arr, FFTW_RODFT11, FFTW_ESTIMATE);
    sphr->plans[1] = FFTWPFX(plan_r2r_1d)(npt, sphr->arr, sphr->arr, FFTW_REDFT11, FFTW_ESTIMATE);
  } else {
    sphr->plans[0] = FFTWPFX(plan_r2r_1d)(npt - 1, sphr->arr+1, sphr->arr+1, FFTW_RODFT00, FFTW_ESTIMATE);
    sphr->plans[1] = FFTWPFX(plan_r2r_1d)(npt + 1, sphr->arr,   sphr->arr,   FFTW_REDFT00, FFTW_ESTIMATE);
  }
#endif

  return sphr;
}

#endif /* defined(FFT) */





#ifdef DHT
/* compute
 *  out(k) = fac Int {from 0 to infinity} dr
 *           r^(D/2) in(r) J_{D/2-1}(k r) */
INLINE void sphr_dht(xdouble *in, xdouble *out, xdouble fac,
    xdht *dht, xdouble *arr, xdouble *r2p, xdouble *k2p)
{
  int i, npt = dht->size;

  for ( i = 0; i < npt; i++ ) arr[i] = in[i] * r2p[i];
  XDHT(apply)(dht, arr, out);
  for ( i = 0; i < npt; i++ ) out[i] *= fac / k2p[i];
}



/* convenience macros */
#define sphr_r2k(sphr, in, out) \
  sphr_dht(in, out, sphr->facr2k, sphr->dht, sphr->arr, sphr->r2p, sphr->k2p)

#define sphr_k2r(sphr, in, out) \
  sphr_dht(in, out, sphr->fack2r, sphr->dht, sphr->arr, sphr->k2p, sphr->r2p)



/* adjust rmax such that r = 1 lies at the middle of the dm'th and dm+1'th bins */
__inline static void sphr_jadjustrmax(sphr_t *sphr, xdouble rmax0)
{
  xdouble dr, km, kp, kM;
  double nu = sphr->dim*.5 - 1;

  dr = rmax0 / sphr->npt;
  sphr->dm = (int)(1/dr + .5);
  kM = gsl_sf_bessel_zero_Jnu(nu, sphr->npt + 1);
  while ( 1 ) {
    km = gsl_sf_bessel_zero_Jnu(nu, sphr->dm);
    kp = gsl_sf_bessel_zero_Jnu(nu, sphr->dm + 1);
    /* adjust rmax such that rmax (j_{nu,k} * .5 + j_{nu,k+1} * .5) / j_{nu,M} = 1 */
    sphr->rmax = kM*2/(km + kp);
    //printf("dm %d, rmax %g, k %g, %g, %g\n", dm, (double) rmax, (double) km, (double) kp, (double) kM);
    if ( sphr->rmax >= rmax0 - 1e-8 || sphr->dm == 1 ) break;
    sphr->dm--;
  }
  sphr->dr = sphr->rmax / sphr->npt;
  printf("D %d, %d bins, %d within the hard core (dr %g), rmax %g, k %g - %g\n",
      sphr->dim, sphr->npt, sphr->dm, (double) sphr->dr,
      (double) sphr->rmax, (double) km, (double) kp);
}



#define sphr_open(dim, npt, rmax, Rmax, ffttype) \
  sphr_opendht(dim, npt, rmax, Rmax, dhtdisk)

__inline static sphr_t *sphr_opendht(int dim, int npt, xdouble rmax,
    xdouble Rmax, int dhtdisk)
{
  int i;
  xdouble tmp1, tmp2, k2;
  sphr_t *sphr;

  sphr = sphr_init(dim, npt);

  if ( FABS(Rmax) > 1e-20 )
    sphr->rmax = Rmax;
  else
    sphr_jadjustrmax(sphr, rmax);

  MAKE1DARR(sphr->r2p,  npt);
  MAKE1DARR(sphr->k2p,  npt);
  sphr->dht = XDHT(newx)(npt, (xdouble) dim/2 - 1,
      sphr->rmax, dhtdisk);
  sphr->surfk *= pow_si(sphr->dht->kmax / sphr->dht->xmax, 2);

  /* find the hard-core boundary dm */
  for ( tmp1 = 0, i = 0; i < npt; i++, tmp1 = tmp2 )
    if ((tmp2 = XDHT(x_sample)(sphr->dht, i)) > 1) {
      printf("r %g, %g, dm %d. \n", (double) tmp1, (double) tmp2, i);
      break;
    }
  sphr->dm = i;

  sphr->facr2k = POW(PI*2, (xdouble) dim/2);
  /* the factor (kmax/xmax)^2 is adapted for the inverse
   * discrete hankel transform using gsl_dht */
  sphr->fack2r = pow_si(sphr->dht->kmax / sphr->dht->xmax, 2);
  sphr->fack2r /= POW(PI*2, (xdouble) dim/2);

  for ( i = 0; i < npt; i++ ) {
    sphr->ri[i]   = XDHT(x_sample)(sphr->dht, i);
    sphr->r2p[i]  = POW(sphr->ri[i], (xdouble) dim/2 - 1);
    sphr->ki[i]   = XDHT(k_sample)(sphr->dht, i);
    sphr->k2p[i]  = POW(sphr->ki[i], (xdouble) dim/2 - 1);

    /* compute r^(dim - 1) dr, used for integration
     *   \int r dr / xmax^2
     * ==>
     *   2/(j_{nu,M})^2 * 1/[J_{nu+1}(j_{nu,k})]^2
     * r = j_{nu,k} / j_{nu,M} */
    k2 = pow_si(sphr->dht->kmax, 2);
    tmp1 = pow_si(sphr->ri[i], dim - 2) * sphr->surfr;
    sphr->rDm1[i] = tmp1 * 2 / (k2 * sphr->dht->J2[i+1]);
    tmp1 = pow_si(sphr->ki[i], dim - 2) * sphr->surfk;
    sphr->kDm1[i] = tmp1 * 2 / (k2 * sphr->dht->J2[i+1]);
  }

  return sphr;
}


#endif /* defined(DHT) */



#endif /* FFTX_H__ */
