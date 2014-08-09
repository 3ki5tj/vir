/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * from the Kirkwood integral equations
 * Test compile (no external library)
 *  gcc ybgvir.c -lm
 * For the normal double precision
 *  gcc -DFFTW -DDHT ybgvir.c -lfftw3 -lgsl -lgslcblas
 * Or for the long double precision
 *  gcc -DFFTW -DDHT -DLDBL ybgvir.c -lfftw3l -lgsl -lgslcblas
 * Or for the 128-bit precision
 *  gcc -DFFTW -DDHT -DF128 ybgvir.c -lfftw3q -lgsl -lgslcblas -lquadmath -lm
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#include "xdouble.h"

#if !defined(FFT0) && !defined(FFTW) && !defined(DHT)
#define FFT0 /* default branch, no external library */
#endif

#ifdef FFT0  /* home-made FFT */
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#endif

#ifdef FFTW /* for odd dimensions */
#include <fftw3.h>
#endif

#if defined(FFT0) || defined(FFTW)
#define FFT
#endif

#ifdef DHT /* for even and odd dimensions */
#include "slowdht.h"
int dhtdisk = SLOWDHT_USEDISK;
#endif


#include "ieutil.h"



#ifndef D
#ifndef DHT
#define D 3
#else
#define D 2
#endif
#endif

int dim = D;
int K; /* (dim - 1) / 2 for an odd dimension */
int nmax = 12;
double rmax = 0;
xdouble Rmax = 0;
int numpt = 32768;
int ffttype = 1;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;
int snapshot = 0;

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
char systitle[32];




#ifdef DHT
#endif /* defined(DHT) */



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the Yvon-Born-Green closure for hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (exact)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type, 0: integral grid points, 1: half-integer grid points");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
#ifdef DHT
  argopt_add(ao, "--disk", "%d", &dhtdisk, "use disk for discrete Hankel transform (DHT)");
  argopt_add(ao, "--tmpdir", "%s", &slowdht_tmpdir, "directory for temporary files");
  argopt_add(ao, "--dhtblock", "%u", &slowdht_block, "block size for DHT input/output");
  argopt_add(ao, "--dhtquiet", "%b", &slowdht_quiet, "suppress DHT output");
#endif
  argopt_add(ao, "-s", "%b", &snapshot, "save intermediate snapshots");
  argopt_add(ao, "-G", "%b", &gaussf, "Gaussian model instead of hard spheres");
  argopt_add(ao, "-I", "%d", &invexp, "exponent of the inverse potential r^(-n)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);

  if ( dim < 1 ) argopt_help(ao);
  K = (dim - 1)/2;
  if ( dim % 2 == 0 ) {
#ifndef DHT
    fprintf(stderr, "cannot do even dimensions %d, define DHT\n", dim);
    exit(1);
#endif
    if ( !argopt_isset(ao, numpt) ) /* adjust the default number of points */
      numpt = 1024;
  }

  if ( rmax <= 0 ) rmax = nmax + 2;
#ifdef FFT
  if ( Rmax <= 0) Rmax = rmax; /* we will adjust it later */
#endif
#ifdef DHT
  if ( Rmax <= 0 ) Rmax = jadjustrmax(rmax, numpt);
#endif

  /* models */
  if ( gaussf || invexp > 0 ) smoothpot = 1;

  { /* write the name of the system */
    char syst[16] = "";
    if ( gaussf ) strcpy(syst, "GF");
    else if ( invexp ) sprintf(syst, "INV%d", invexp);
    sprintf(systitle, "%s%s", (dim % 2 ? "" : "h"), syst);
  }

  printf("D %d, Rmax %f, npt %d, %s\n",
      dim, (double) Rmax, numpt, systitle);

  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



__inline static xdouble *summ(int l, int npt, xdouble **in, xdouble *out)
{
  int i, m;

  for ( i = 0; i < npt; i++ )
    for ( out[i] = 0, m = 0; m <= l; m++ )
      out[i] += in[m][i];
  return out;
}



/* compute w(k) from c(k) and h(k) */
__inline static void get_wk_kirk(int l, int npt, xdouble **wkl,
    xdouble ***ck, xdouble **hk)
{
  int i, u, m;

  for ( m = 1; m <= l; m++ ) { /* ksi^m */
    for ( i = 0; i < npt; i++ )
      wkl[m][i] = 0;
    /* w(k) = rho c(k) h(k)
     * w_{l, m}(k) = Sum_{u from 0 to l - 1} c_{u, m}(k) h_{l - u - 1}(k) */
    for ( u = 0; u < l; u++ )
      for ( i = 0; i < npt; i++ )
        wkl[m][i] += ck[u][m-1][i] * hk[l - u - 1][i];
  }
}



/* y(r) = exp(w(r)) */
__inline static void get_yr_kirk(int l, int npt, xdouble ***yr, xdouble ***wr)
{
  int i, m, l1, m1, m1b;
  xdouble y;

  for ( i = 0; i < npt; i++ ) {
    for ( m = 1; m <= l; m++ ) {
      y = 0;
      for ( l1 = 1; l1 <= l; l1++ ) {
        /* m1 >= max{1, m - l + l1}, l - l1 >= m - m1 */
        if ( (m1 = m - l + l1) < 1 ) m1 = 1;
        /* m1 <= min{m, l1} */
        if ( (m1b = m) > l1 ) m1b = l1;
        /* dy/d(rho) = dw/d(rho) y
         * Sum_{l, m} l y_{l, m} rho^{l-1} ksi^m
         * = Sum_{l1, m1} l1 w_{l1, m1} rho^{l1-1} ksi^{m1}
         *   Sum_{l2, m2}    w_{l2, m2} rho^{l2}   ksi^{m2} */
        for ( m1 = 1; m1 <= m1b; m1++ )
          y += l1 * wr[l1][m1][i] * yr[l-l1][m-m1][i];
      }
      yr[l][m][i] = y / l;
    }
  }
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype)
{
  xdouble facr2k, fack2r, surfr, surfk;
  xdouble Bc, Bv, By = 0, B2, fcorr = 0;
  xdouble *fr, *fk, *rdfr = NULL;
  xdouble ***wr = NULL, *hrl = NULL, *crl = NULL, ***yr = NULL;
  xdouble **wkl = NULL, ***ck = NULL, **hk = NULL, *hk0;
  xdouble *arr;
  xdouble *ri, *ki, *rDm1, *kDm1;
  int i, dm, l, l0 = 1, m;
  clock_t t0 = clock(), t1;


#ifdef FFT /* for odd dimensions */
  xdouble dr,dk, rl, invrl, kl, invkl;
  xdouble **r2pow = NULL, **invr2pow = NULL, **k2pow = NULL, **invk2pow = NULL;
  long *coef = NULL;
  FFTWPFX(plan) plans[2] = {NULL, NULL};

  rmax = adjustrmax(rmax, numpt, &dr, &dm, ffttype);
  dk = PI/dr/npt;
  facr2k = pow_si(PI*2, K) *  dr;
  fack2r = pow_si(PI*2, -K-1) * dk;
#endif /* defined(FFT) */


#ifdef DHT /* for even dimensions */
  xdouble tmp1, tmp2;
  xdouble *r2p, *k2p;
  xdht *dht = NULL;

  dht = XDHT(newx)(npt, (xdouble) dim/2 - 1, rmax, dhtdisk);

  for ( tmp1 = 0, dm = 0; dm < npt; dm++, tmp1 = tmp2 )
    if ((tmp2 = XDHT(x_sample)(dht, dm)) > 1) {
      printf("r %g, %g, dm %d. \n", (double) tmp1, (double) tmp2, dm);
      break;
    }

  facr2k = POW(PI*2, (xdouble) dim/2);
  /* the factor (kmax/xmax)^2 is adapted for the inverse
   * discrete hankel transform using gsl_dht */
  fack2r = pow_si(dht->kmax / dht->xmax, 2);
  fack2r /= POW(PI*2, (xdouble) dim/2);
#endif /* defined(DHT) */


  /* B2 = (1/2) (PI*2)^(dim/2) / dim!! for an even dim
   *    = (PI*2)^((dim - 1)/2) / dim!! for an odd dim */
  B2 = dim % 2 ? 1 : 0.5;
  for ( i = 2 + dim % 2; i <= dim; i += 2 ) B2 *= PI*2/i;
  surfr = B2 * 2 * dim;
  surfk = surfr / pow_si(PI*2, dim);
#ifdef DHT
  surfk *= pow_si(dht->kmax / dht->xmax, 2);
#endif /* defined(DHT) */
  if ( gaussf ) B2 = SQRT(pow_si(PI, dim))/2;
  /* TODO: compute B2 for inverse potential */
  if ( invexp > 0 ) B2 = 0;
  printf("D %d, dr %f, dm %d, rmax %f, ffttype %d, B2 %.10e\n",
      dim, (double) rmax/npt, dm, (double) rmax, ffttype, (double) B2);

  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE1DARR(rDm1, npt); /* r^(D - 1) dr */
  MAKE1DARR(kDm1, npt); /* k^(D - 1) dk */


#ifdef FFT
  if ( K > 0 ) {
    MAKE2DARR(r2pow, K, npt) /* r^{K - l} for l = 0, ..., K */
    MAKE2DARR(k2pow, K, npt) /* k^{K - l} for l = 0, ..., K */
    MAKE2DARR(invr2pow, K, npt) /* r^{-K-l} for l = 0, ..., K */
    MAKE2DARR(invk2pow, K, npt) /* k^{-K-l} for l = 0, ..., K */
  }
  for ( i = 0; i < npt; i++ ) {
    ri[i]  = dr * (i*2 + (ffttype ? 1 : 0))/2;
    ki[i]  = dk * (i*2 + (ffttype ? 1 : 0))/2;
    rl = 1;
    kl = 1;
    for ( l = 1; l <= K; l++ ) {
      r2pow[K - l][i] = (rl *= ri[i]);
      k2pow[K - l][i] = (kl *= ki[i]);
    }

    if ( ffttype == 0 && i == 0 ) continue;

    invrl = 1/rl;
    invkl = 1/kl;
    for ( l = 0; l < K; l++ ) {
      invr2pow[l][i] = invrl;
      invk2pow[l][i] = invkl;
      invrl /= ri[i];
      invkl /= ki[i];
    }

    rDm1[i] = surfr * pow_si(ri[i], dim - 1) * dr;
    kDm1[i] = surfk * pow_si(ki[i], dim - 1) * dk;
  }
#endif /* defined(FFT) */


#ifdef DHT
  MAKE1DARR(r2p,  npt);
  MAKE1DARR(k2p,  npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i]   = XDHT(x_sample)(dht, i);
    r2p[i]  = POW(ri[i], (xdouble) dim/2 - 1);
    ki[i]   = XDHT(k_sample)(dht, i);
    k2p[i]  = POW(ki[i], (xdouble) dim/2 - 1);

    /* compute r^(dim - 1) dr, used for integration
     *   \int r dr / xmax^2
     * ==>
     *   2/(j_{nu,M})^2 * 1/[J_{nu+1}(j_{nu,k})]^2
     * r = j_{nu,k} / j_{nu,M} */
    tmp1 = pow_si(ri[i], dim - 2) * surfr;
    rDm1[i] = tmp1 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
    tmp1 = pow_si(ki[i], dim - 2) * surfk;
    kDm1[i] = tmp1 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
  }
#endif /* defined(DHT) */


  /* auxiliary array for FFTW or DHT
   * needs npt + 1 elements for FFTW_REDFT00 */
  MAKE1DARR(arr, npt + 1);

#ifdef FFT
  /* compute the coefficients of the spherical Bessel function */
  if ( K > 0 ) {
    MAKE1DARR(coef, K);
    getjn(coef, K - 1);
  }
#endif

#ifdef FFTW
  /* plans[0] is the sine transform, plans[1] is the cosine transform */
  if ( ffttype ) {
    plans[0] = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);
    plans[1] = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_REDFT11, FFTW_ESTIMATE);
  } else {
    plans[0] = FFTWPFX(plan_r2r_1d)(npt - 1, arr + 1, arr + 1, FFTW_RODFT00, FFTW_ESTIMATE);
    plans[1] = FFTWPFX(plan_r2r_1d)(npt + 1, arr, arr, FFTW_REDFT00, FFTW_ESTIMATE);
  }
#endif

  MAKE1DARR(fr, npt);
  MAKE1DARR(fk, npt);
  if ( smoothpot ) MAKE1DARR(rdfr, npt);
  MAKE1DARR(hrl, npt);
  MAKE1DARR(crl, npt);
  MAKE2DARR(wkl, nmax - 1, npt);
  MAKE2DARR(hk,  nmax - 1, npt);
  MAKE1DARR(hk0, nmax - 1);
  MAKE3DARR(wr,  nmax - 1, nmax - 1, npt);
  MAKE3DARR(yr,  nmax - 1, nmax - 1, npt);
  MAKE3DARR(ck,  nmax - 1, nmax - 1, npt);

  hk0[0] = -2*B2;

  /* compute f(r) */
  for ( i = 0; i < npt; i++ ) { /* compute f(r) = exp(-beta u(r)) - 1 */
    if ( gaussf ) { /* f(r) = exp(-r^2) */
      xdouble r2 = ri[i] * ri[i];
      fr[i] = -EXP(-r2);
      rdfr[i] = -2*r2*fr[i];
    } else if ( invexp > 0 ) { /* inverse potential r^(-invexp) */
      xdouble pot = POW(ri[i], -invexp);
      fr[i] = EXP(-pot) - 1;
      rdfr[i] = pot * invexp * (fr[i] + 1);
    } else { /* hard-sphere */
      fr[i] = (i < dm) ? -1 : 0;
    }
    hrl[i] = fr[i];
  }
  /* compute fk */
  sphr_r2k(fr, fk);
  COPY1DARR(hk[0], fk, npt); /* hk[0] = fk */
  COPY1DARR(ck[0][0], fk, npt);

  for ( i = 0; i < npt; i++ )
    yr[0][0][i] = 1;

  t1 = clock();

  fnvir = savevirheadx(fnvir, systitle, dim, l0, nmax,
      IETYPE_KIRKWOOD, 0, 0, npt, rmax, t1 - t0, 1, 1, -1, 0, 0, 0);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* get w_{l,m}(k) */
    get_wk_kirk(l, npt, wkl, ck, hk);

    /* w_{l,m}(k) --> w_{l,m}(r) */
    for ( m = 1; m <= l; m++ )
      sphr_k2r(wkl[m], wr[l][m]);

#ifdef PY
    for ( m = 1; m <= l; m++ )
      for ( i = 0; i < npt; i++ )
        yr[l][m][i] = wr[l][m][i];
#else
    /* y_{l,m}(r) = [exp w(r)]_{l,m} */
    get_yr_kirk(l, npt, yr, wr);
#endif

    for ( i = 0; i < npt; i++ ) {
      /* h(r) = (f(r) + 1) y(r) - 1
       * h_l(r) = (f(r) + 1) Sum_{m = 1 to l} y_{l,m}(r) (with ksi = 1)  */
      hrl[i] = 0;
      for ( m = 1; m <= l; m++ )
        hrl[i] += (fr[i] + 1) * yr[l][m][i];
    }

    /* h_l(r) --> h_l(k) */
    sphr_r2k(hrl, hk[l]);

    /* compute c_{l, m+1}(r) = y_{l, m}(r) f(r) / (m+1) */
    for ( m = 1; m <= l; m++ ) {
#ifdef HNC
      /* this is the HNC case */
      for ( i = 0; i < npt; i++)
        crl[i] = (fr[i] + 1) * yr[l][m][i] - wr[l][m][i];
#elif defined(PY)
      /* this is the PY case */
      for ( i = 0; i < npt; i++)
        crl[i] = wr[l][m][i] * fr[i];
#else
      for ( i = 0; i < npt; i++ )
        crl[i] = yr[l][m][i] * fr[i] / (m + 1);
      /* c_{l, m+1}(r) --> c_{l, m+1}(k) */
#endif
      sphr_r2k(crl, ck[l][m]);
    }

    summ(l, npt, yr[l], crl);
    Bv = get_Bv(npt, crl, smoothpot, rdfr, rDm1, dim, dm, B2);
    Bc = get_Bc_hr(l, npt, hrl, hk0, rDm1);
    summ(l, npt, wr[l], crl);
    By = get_zerosep(crl, ri)*l/(l+1);

    savevir(fnvir, dim, l+2, Bc, Bv, 0, 0, 0, By, B2, 0, fcorr);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(arr,  npt);
  FREE1DARR(ri,   npt);
  FREE1DARR(ki,   npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr,   npt);
  FREE1DARR(fk,   npt);
  FREE1DARR(rdfr, npt);

  FREE1DARR(hrl, npt);
  FREE1DARR(crl, npt);
  FREE2DARR(wkl, nmax - 1, npt);
  FREE2DARR(hk,  nmax - 1, npt);
  FREE1DARR(hk0, nmax - 1);
  FREE3DARR(wr,  nmax - 1, nmax - 1, npt);
  FREE3DARR(yr,  nmax - 1, nmax - 1, npt);
  FREE3DARR(ck,  nmax - 1, nmax - 1, npt);

#ifdef FFT
  FREE1DARR(coef, K);
#ifdef FFTW
  FFTWPFX(destroy_plan)(plans[0]);
  FFTWPFX(destroy_plan)(plans[1]);
#endif /* defined FFTW */
  FREE2DARR(r2pow, K, npt);
  FREE2DARR(k2pow, K, npt);
  FREE2DARR(invr2pow, K, npt);
  FREE2DARR(invk2pow, K, npt);
#endif /* defined(FFT) */

#ifdef DHT
  FREE1DARR(r2p,  npt);
  FREE1DARR(k2p,  npt);
  XDHT(free)(dht);
#endif /* defined(DHT) */


  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, Rmax, ffttype);
#ifdef FFTW
  FFTWPFX(cleanup)();
#endif
  return 0;
}
