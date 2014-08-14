/* Computing the virial coefficients of an hard-sphere fluid
 * from the PY, HNC and other integral equations
 * Test compile (no external library)
 *  gcc ievir.c -lm
 * For the normal double precision
 *  gcc -DFFTW -DDHT ievir.c -lfftw3 -lgsl -lgslcblas
 * Or for the long double precision
 *  gcc -DFFTW -DDHT -DLDBL ievir.c -lfftw3l -lgsl -lgslcblas
 * Or for the 128-bit precision
 *  gcc -DFFTW -DDHT -DF128 ievir.c -lfftw3q -lgsl -lgslcblas -lquadmath -lm
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
int dohnc = 0;
int singer = 0;
int ring = 0;
int mkcorr = 0;
int expcorr = 0; /* exponential correction */
int fast = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;
int snapshot = 0;

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
char systitle[32];

int ietype = 0;

xdouble hncamp = 0, hncq = 1; /* Marucho-Pettitt */
xdouble hncalpha = -1; /* Rogers-Young */
xdouble hcs = 1; /* Hutchinson-Conkie s */
xdouble rowphi = 0; /* Rowlinson's Phi */
xdouble invphi = 0; /* inverse Rowlinson's Phi */
xdouble hurstm = 0; /* Hurst's m */
xdouble verleta = 0, verletb = 0; /* Verlet modified */
xdouble bbpgs = 15./8; /* MS/BBPG s */
xdouble sqrs = 0.5;
xdouble xsqrs = 0;
xdouble cubs = 0.5;
xdouble cubt = 0.0;
xdouble geos = 1;
xdouble logs = 1;

xdouble shift = 0, shiftinc = 0;
int shiftl0 = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  xdouble shiftn = 0;
  int usems = 0;

  ao->desc = "computing the virial coefficients from the PY/HNC closure for hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (exact)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type, 0: integral grid points, 1: half-integer grid points");
  argopt_add(ao, "-T", "%d", &ietype, "type of closure, 0: PY, 1: HNC-like; set the respective parameters automatically sets the type");
  argopt_add(ao, "-q", "%" XDBLSCNF "f", &hncq,  "q of the hypernetted-chain (Marucho-Pettitt) approximation, y(r) = a0 exp(q t(r)) + (1-a0) + (1-a0 q) t(r)");
  argopt_add(ao, "-a", "%" XDBLSCNF "f", &hncamp, "a0 of the hypernetted-chain (Marucho-Pettitt) approximation, 0: to be set by 1/q");
  argopt_add(ao, "--rya", "%" XDBLSCNF "f", &hncalpha, "alpha, in the Roger-Young switch function [1 - exp(-alpha*r)]^m");
  argopt_add(ao, "--hcs", "%" XDBLSCNF "f", &hcs, "s, in the Hutchinson-Conkie closure, y(r) = (1 + s t(r))^(1/s)");
  argopt_add(ao, "--rphi", "%" XDBLSCNF "f", &rowphi, "phi of the Rowlinson approximation, t(r) = (1 - phi) (y(r) - 1) + phi ln y(r)");
  argopt_add(ao, "--iphi", "%" XDBLSCNF "f", &invphi, "phi of the inverse Rowlinson approximation, y(r) = exp(t(r)) phi + (1 - phi) (1 + t(r))");
  argopt_add(ao, "--hm", "%" XDBLSCNF "f", &hurstm, "m of the Hurst approximation, t(r) = [y(r)]^m log y(r)");
  argopt_add(ao, "--va", "%" XDBLSCNF "f", &verleta, "A of the Verlet approximation, y(r) = exp[t(r) - A t(r)^2 /2 / (1 + B t(r) / 2) ]");
  argopt_add(ao, "--vb", "%" XDBLSCNF "f", &verletb, "B of the Verlet approximation");
  argopt_add(ao, "--bbpgs", "%" XDBLSCNF "f", &bbpgs, "s of the BBPG approximation, y(r) = exp[ (1 + s t(r))^(1/s) - 1 ]");
  argopt_add(ao, "--ms", "%b", &usems, "Martynov-Sarkisov approximation, y(r) = exp[ (1 + 2 t(r))^(1/2) - 1 ]");
  argopt_add(ao, "--sqrs", "%" XDBLSCNF "f", &sqrs, "s of the quadratic approximation, y(r) = 1 + t(r) + s t(r)^2 / 2");
  argopt_add(ao, "--xsqrs", "%" XDBLSCNF "f", &xsqrs, "s of the exponential quadratic approximation, y(r) = exp[ t(r) + s t(r)^2 ]");
  argopt_add(ao, "--cubs", "%" XDBLSCNF "f", &cubs, "s of the cubic approximation, y(r) = 1 + t(r) + s t(r)^2 / 2 + t t(r)^3/6");
  argopt_add(ao, "--cubt", "%" XDBLSCNF "f", &cubt, "t of the cubic approximation");
  argopt_add(ao, "--geos", "%" XDBLSCNF "f", &geos, "s of the geometric approximation, y(r) = 1 + t(r) + (s/2) t(r)^2 + (s/2)^2 t(r)^3 + ...");
  argopt_add(ao, "--logs", "%" XDBLSCNF "f", &logs, "s of the logarithmic approximation, y(r) = 1 + t(r) + s t(r)^2/2 + s^2 t(r)^3/3 + ...");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted-chain (HNC) like approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "linearly correct the closure");
  argopt_add(ao, "--expc", "%b", &expcorr, "exponentially correct the closure");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--fast", "%b", &fast, "save tk and yr to accelerate the calculation");
  argopt_add(ao, "-c", "%" XDBLSCNF "f", &shift, "shift of t(r) in computing Bv");
  argopt_add(ao, "-C", "%" XDBLSCNF "f", &shiftn, "shift of t(r) in computing Bv, in terms of n");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &shiftinc, "increment of the shift");
  argopt_add(ao, "-L", "%d", &shiftl0, "minimal order l for the shift");
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
  if ( Rmax <= 0 ) Rmax = jadjustrmax(rmax, numpt, dim);
#endif

  /* closure type */
  if ( hncq > 0 && hncamp <= 0 ) /* Marucho-Pettitt extended HNC */
    hncamp = 1/hncq; /* set amp automatically = 1/q */
  if ( argopt_isset(ao, hncamp) || argopt_isset(ao, hncq)
    || argopt_isset(ao, hncalpha) ) {
    ietype = IETYPE_HNC;
    dohnc = 1; /* turn on HNC, if necessary */
  }

  if ( argopt_isset(ao, hcs) ) ietype = IETYPE_HC;
  if ( argopt_isset(ao, rowphi) ) ietype = IETYPE_ROWLINSON;
  if ( argopt_isset(ao, invphi) ) ietype = IETYPE_INVROWLINSON;
  if ( argopt_isset(ao, hurstm) ) ietype = IETYPE_HURST;
  if ( argopt_isset(ao, verleta) || argopt_isset(ao, verletb) ) ietype = IETYPE_VERLET;
  if ( argopt_isset(ao, usems) ) { ietype = IETYPE_BBPG; bbpgs = 2; }
  if ( argopt_isset(ao, bbpgs) ) ietype = IETYPE_BBPG;
  if ( argopt_isset(ao, sqrs) ) ietype = IETYPE_SQR;
  if ( argopt_isset(ao, xsqrs) ) ietype = IETYPE_XSQR;
  if ( argopt_isset(ao, cubs) || argopt_isset(ao, cubt) ) ietype = IETYPE_CUB;
  if ( argopt_isset(ao, geos) ) ietype = IETYPE_GEO;
  if ( argopt_isset(ao, logs) ) ietype = IETYPE_LOG;

  if ( ietype > IETYPE_HNC ) dohnc = 0;

  if ( dohnc ) ietype = IETYPE_HNC;

  if ( argopt_isset(ao, shift) || argopt_isset(ao, shiftn) || argopt_isset(ao, shiftinc) )
    mkcorr = 1;

  if ( mkcorr || ietype > IETYPE_HNC ) /* Singer and ring formulas are inapplicable to corrections */
    singer = ring = 0;
  if ( singer ) ring = 1;

  /* correction parameters */
  if ( mkcorr && ietype == IETYPE_PY ) { dohnc = 1; }
  if ( argopt_isset(ao, shiftn) ) shift = shiftn + shiftinc*2;
  /* for the Gaussian model, Bv is exact for l <= 2 (n <= 4)
   * otherwise, Bv is exact for l = 1 (n = 3) */
  if ( shiftl0 <= 0 ) shiftl0 = gaussf ? 3 : 2;

  /* models */
  if ( gaussf || invexp > 0 ) smoothpot = 1;

  { /* write the name of the system */
    char syst[16] = "";
    if ( gaussf ) strcpy(syst, "GF");
    else if ( invexp ) sprintf(syst, "INV%d", invexp);
    sprintf(systitle, "%s%s", (dim % 2 ? "" : "h"), syst);
  }

  printf("D %d, Rmax %f, npt %d, ietype %d, HNC %d, a0 %g, q %g, %s\n",
      dim, (double) Rmax, numpt, ietype, dohnc,
      (double) hncamp, (double) hncq, systitle);

  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype)
{
  xdouble facr2k, fack2r, surfr, surfk;
  xdouble Bc, Bv, Bm = 0, Bh = 0, Br = 0, By = 0, B2, fcorr = 0;
  xdouble *fr, *rdfr = NULL, *swr = NULL;
  xdouble *crl, *trl, *yrl, *tkl;
  xdouble **ck, **tk = NULL, **cr = NULL, **tr = NULL, **yr = NULL;
  xdouble *arr, *yrcoef = NULL;
  xdouble *vc = NULL, *brl = NULL, **br = NULL; /* correction functions */
  xdouble *ri, *ki, *rDm1, *kDm1;
  xdouble *yr0, *lnyr0;
  int i, dm, l, l0 = 1;
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
  printf("D %d, dr %f, dm %d, rmax %f, ffttype %d, ietype %d, B2 %.10e\n",
      dim, (double) rmax/npt, dm, (double) rmax, ffttype, ietype, (double) B2);

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
  if ( smoothpot ) MAKE1DARR(rdfr, npt);
  MAKE1DARR(swr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE1DARR(yrl, npt);
  MAKE1DARR(tkl, npt);
  MAKE2DARR(ck, nmax - 1, npt);
  if ( fast ) {
    MAKE2DARR(tk, nmax - 1, npt);
  }

  /* compute f(r) and f(k) = c0(k) */
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
    crl[i] = fr[i];
  }

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt);
  }

  if ( ietype != IETYPE_PY || mkcorr || expcorr || singer ) {
    MAKE2DARR(tr, nmax - 1, npt);
  }

  MAKE1DARR(yr0, nmax - 1);
  MAKE1DARR(lnyr0, nmax - 1);
  yr0[0] = 1;
  lnyr0[0] = 0;

  if ( ietype > IETYPE_HNC ) {
    /* compute the coefficients, they are not actually used
     * but they offer an alternative (series expansion) as a backup */
    MAKE1DARR(yrcoef, nmax - 1);
    if ( ietype == IETYPE_HC ) {
      init_hccoef(yrcoef, nmax - 1, hcs);
    } else if ( ietype == IETYPE_BBPG ) {
      init_bbpgcoef(yrcoef, nmax - 1, bbpgs);
    } else if ( ietype == IETYPE_ROWLINSON ) {
      init_rowlinsoncoef(yrcoef, nmax - 1, rowphi);
    } else if ( ietype == IETYPE_INVROWLINSON ) {
      init_invrowlinsoncoef(yrcoef, nmax - 1, invphi);
    } else if ( ietype == IETYPE_HURST ) {
      init_hurstcoef(yrcoef, nmax - 1, hurstm);
    } else if ( ietype == IETYPE_VERLET ) {
      init_verletcoef(yrcoef, nmax - 1, verleta, verletb);
    } else if ( ietype == IETYPE_SQR ) {
      init_sqrcoef(yrcoef, nmax - 1, sqrs);
    } else if ( ietype == IETYPE_XSQR ) {
      init_xsqrcoef(yrcoef, nmax - 1, xsqrs);
    } else if ( ietype == IETYPE_CUB ) {
      init_cubcoef(yrcoef, nmax - 1, cubs, cubt);
    } else if ( ietype == IETYPE_GEO ) {
      init_geocoef(yrcoef, nmax - 1, geos);
    } else if ( ietype == IETYPE_LOG ) {
      init_logcoef(yrcoef, nmax - 1, logs);
    }
    print_yrcoef(yrcoef, 7);
  }

  if ( dohnc ) { /* initialize the Rogers-Young switch function */
    for ( i = 0; i < npt; i++ ) {
      swr[i] = 1;
      if ( hncalpha >= 0 )
        swr[i] = 1 - EXP( -hncalpha * ri[i] );
    }
  }

  if ( fast ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) yr[0][i] = 1;
  }

  if ( mkcorr ) {
    MAKE1DARR(brl, npt);
  }

  if ( expcorr ) { /* allocate space for the bridge function */
    MAKE2DARR(br, nmax - 1, npt);
  }

  t1 = clock();

  if ( snapshot ) /* TODO: snapshot code is obsolete */
    l0 = snapshot_open(dim, nmax, rmax, dohnc, mkcorr, ring, singer,
        npt, ck, tk, cr, tr, crl, trl, yr);

  fnvir = savevirheadx(fnvir, systitle, dim, l0, nmax,
      ietype, mkcorr, expcorr, npt, rmax, t1 - t0,
      hncamp, hncq, hncalpha, shift, shiftinc, shiftl0);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* c_l(r) --> c_l(k) for the previous l */
    sphr_r2k(crl, ck[l-1]);

    /* compute the ring sum based on ck */
    if ( ring ) {
      Bh = get_BhBrk(l, npt, dohnc, ck, kDm1, ki, &Br);
    }

    /* compute t_l(k) from c_0(k), ... c_{l-1}(k) */
    if ( tk != NULL ) { /* fast version requires tk[] */
      get_tk_oz_fast(l, npt, ck, tk);
      COPY1DARR(tkl, tk[l], npt); /* tkl = tk[l] */
    } else { /* slow version requires ck[] only */
      get_tk_oz(l, npt, ck, tkl);
    }

    /* t_l(k) --> t_l(r) */
    sphr_k2r(tkl, trl);

    if ( tr != NULL ) { /* save tr if needed */
      COPY1DARR(tr[l], trl, npt);
    }

    /* compute the cavity function y(r) */
    if ( expcorr ) {
      get_yr_hnc_br(l, npt, yrl, tr, br);
    } else
    if ( ietype == IETYPE_PY && !mkcorr ) { /* PY closure */
      for ( i = 0; i < npt; i++ )
        yrl[i] = trl[i];
    } else if ( ietype == IETYPE_HNC
            || (ietype == IETYPE_PY && mkcorr) ) { /* HNC closure */
      if ( yr != NULL ) {
        get_yr_hncx_fast(l, npt, yrl, yr, tr, hncamp, hncq, swr);
      } else {
        get_yr_hncx(l, npt, yrl, tr, hncamp, hncq, swr);
      }
    } else if ( ietype == IETYPE_HC ) {
      get_yr_hc(l, npt, yrl, tr, hcs);
    } else if ( ietype == IETYPE_BBPG ) {
      get_yr_bbpg(l, npt, yrl, tr, bbpgs);
    } else if ( ietype == IETYPE_INVROWLINSON ) {
      get_yr_invrowlinson(l, npt, yrl, tr, invphi);
    } else if ( ietype == IETYPE_VERLET ) {
      get_yr_verlet(l, npt, yrl, tr, verleta, verletb);
    } else if ( ietype == IETYPE_SQR ) {
      get_yr_sqr(l, npt, yrl, tr, sqrs);
    } else if ( ietype == IETYPE_XSQR ) {
      get_yr_xsqr(l, npt, yrl, tr, xsqrs);
    } else if ( ietype == IETYPE_CUB ) {
      get_yr_cub(l, npt, yrl, tr, cubs, cubt);
    } else if ( ietype == IETYPE_GEO ) {
      get_yr_geo(l, npt, yrl, tr, geos);
    } else if ( ietype == IETYPE_LOG ) {
      get_yr_log(l, npt, yrl, tr, logs);
    } else if ( ietype == IETYPE_ROWLINSON ) {
      get_yr_rowlinson(l, npt, yrl, tr, rowphi);
    } else if ( ietype == IETYPE_HURST ) {
      get_yr_hurst(l, npt, yrl, tr, hurstm);
    } else { /* use it as the last resort or a check */
      get_yr_series(l, npt, yrl, tr, yrcoef);
    }

    yr0[l] = get_zerosep(mkcorr ? trl : yrl, ri);

    for ( i = 0; i < npt; i++ ) {
      /* c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      /* in case we need to make a correction,
       * we still use t(r) as y(r), and y(r) is used in
       * constructing the correction */
      crl[i] = (fr[i] + 1) * (mkcorr ? trl[i] : yrl[i]) - trl[i];
    }

    Bv = get_Bv(npt, mkcorr ? trl : yrl, smoothpot, rdfr, rDm1, dim, dm, B2);
    Bc = get_Bc(l, npt, crl, rDm1);

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( dohnc ) {
        Bm = get_Bm_singer(l, npt, cr, tr, rDm1);
        Bh += get_Bh_singer(l, npt, cr, tr, rDm1);
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, rDm1);
      }
    } else {
      Bm = Bh = 0;
    }

    if ( mkcorr || expcorr ) {
      vc = mkcorr ? brl : br[l];

      if ( mkcorr ) {
        /* construct the linear correction function */
        for ( i = 0; i < npt; i++ )
          vc[i] = yrl[i] - trl[i];
      } else {
        /* construct the exponential correction function */
        get_exp_corr(l, npt, vc, tr);
      }

      /* apply the shift */
      if ( l >= shiftl0 ) Bv *= 1 + shift + shiftinc * l;

      Bm = get_corr1x(l, npt, dm, crl, fr, rdfr, rDm1,
                      dim, B2, vc, &Bc, &Bv, &fcorr);
      yr0[l] += get_zerosep(vc, ri);
    }

    By = update_lnyr0(l, yr0, lnyr0);
    savevir(fnvir, dim, l+2, Bc, Bv, Bm, Bh, Br, By, B2, mkcorr|expcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yrl);
    if ( snapshot )
      snapshot_take(l, npt, ck[l-1], tk[l], crl, trl, nmax, yr);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(arr,  npt);
  FREE1DARR(ri,   npt);
  FREE1DARR(ki,   npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr,   npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(swr,  npt);
  FREE1DARR(crl,  npt);
  FREE1DARR(trl,  npt);
  FREE1DARR(yrl,  npt);
  FREE1DARR(tkl,  npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE2DARR(br, nmax - 1, npt);
  FREE1DARR(brl, npt);
  FREE1DARR(yrcoef, nmax - 1);
  FREE1DARR(yr0, nmax - 1);
  FREE1DARR(lnyr0, nmax - 1);


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
