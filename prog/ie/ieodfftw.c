/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ieodfftw.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ieodfftw.c -lfftw3l
 * Or for the 128-bit precision
 *  gcc -DF128 ieodfftw.c -lfftw3q -lquadmath -lm
 * To disable FFTW
 *  gcc -DNOFFTW ieodfftw.c -lm
 * */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#include "xdouble.h"

#ifdef NOFFTW
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#else
#include <fftw3.h>
#endif

#include "ieutil.h"



#ifdef D
int dim = D;
#else
int dim = 3;
#endif

int K;
int nmax = 12;
double rmax = 0;
int numpt = 32768;
int ffttype = 1;
int dohnc = 0;
int singer = 0;
int ring = 0;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;
int snapshot = 0;

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
char systitle[32];

xdouble hncamp = 1, hncq = 1, hncalpha = -1;
xdouble shift = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &dim, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", "%lf", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted chain approximation");
  argopt_add(ao, "-a", "%" XDBLSCNF "f", &hncamp, "a0 of the hypernetted chain approximation");
  argopt_add(ao, "-q", "%" XDBLSCNF "f", &hncq,  "q  of the hypernetted chain approximation");
  argopt_add(ao, "-A", "%" XDBLSCNF "f", &hncalpha, "alpha, in the Roger-Young switch function [1 - exp(-alpha*r)]^m");
  argopt_add(ao, "-c", "%" XDBLSCNF "f", &shift, "shift of t(r) in computing Bv");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "-s", "%b", &snapshot, "save intermediate snapshots");
  argopt_add(ao, "-G", "%b", &gaussf, "Gaussian model instead of hard spheres");
  argopt_add(ao, "-I", "%d", &invexp, "exponent of the inverse potential r^(-n)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( dim < 3 || dim % 2 == 0 ) argopt_help(ao);
  K = (dim - 1)/2;
  if ( rmax <= 0 ) rmax = nmax + 2;
  if ( mkcorr )
    singer = ring = 0;
  else if ( FABS(hncamp - 1) > 1e-6 || FABS(hncq - 1) > 1e-6 || hncalpha >= 0 )
    dohnc = 1; /* turn on HNC, if necessary */
  if ( FABS(hncq - 1) > 1e-6 && hncq > 0 && hncamp <= 0 )
    hncamp = 1/hncq; /* set amp automatically = 1/q */
  if ( singer ) ring = 1;
  if ( gaussf || invexp > 0 ) smoothpot = 1;
  if ( gaussf ) strcpy(systitle, "GF");
  else if ( invexp ) sprintf(systitle, "INV%d", invexp);
  printf("D %d, K %d, rmax %f, npt %d, HNC %d, %s\n",
      dim, K, (double) rmax, numpt, dohnc, systitle);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}




/* get the coefficient c_l of spherical Bessel function jn(x)
 *  jn(x) = Sum_{l = 0 to n} c_l [sin(x) or cos(x)] / x^{n + l + 1},
 * where the l'th term is sin(x) if n + l is even, or cos(x) otherwise */
static void getjn(long *c, int n)
{
  int i, k;
  const char *fs[2] = {"sin(x)", "cos(x)"};

  c[0] = 1; /* j0 = sin(x)/x; */
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
    printf(" %+ld*%s/x^%d", c[i], fs[(i+n)%2], n + i + 1);
  printf("\n");
}



/* compute
 *    out(k) = 2 fac0 Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' are used in the intermediate steps by the FFTW plans `p'
 * */
static void sphr(int npt, xdouble *in, xdouble *out, xdouble fac,
    FFTWPFX(plan) p[2], xdouble *arr, long *coef,
    xdouble **r2p, xdouble **k2q, int ffttype)
{
  int i, l, iscos;

  if ( ffttype != 0 && ffttype != 1 ) return;

  /* clear the output */
  for ( i = 0; i < npt; i++ ) out[i] = 0;

  /* several rounds of transforms */
  for ( l = 0; l < K; l++ ) {
    /* decide if we want to do a sine or cosine transform */
    iscos = (K + l + 1) % 2;
    for ( i = 0; i < npt; i++ ) {
      /* arr[i] = in[i] * pow(dx*(2*i + 1)/2, K - l) * coef[l]; */
      arr[i] = in[i] * coef[l] * r2p[l][i];
    }
#ifdef NOFFTW
    if ( iscos ) {
      if ( ffttype ) {
        cost11(arr, npt);
      } else {
        cost00(arr, npt);
      }
    } else {
      if ( ffttype ) {
        sint11(arr, npt);
      } else {
        sint00(arr, npt);
      }
    }
#else
    FFTWPFX(execute)(p[iscos]);
#endif
    if ( iscos && ffttype == 0 ) arr[npt] = 0;
    for ( i = 0; i < npt; i++ ) {
      /* out[i] += arr[i] * fac / pow(dk*(2*i + 1)/2, K + l); */
      out[i] += arr[i] * fac * k2q[l][i];
    }
  }
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype, int dohnc)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk;
  xdouble Bc, Bv, Bm = 0, Bh = 0, Br = 0, B2, fcorr = 0;
  xdouble *fr, *rdfr = NULL, *swr = NULL;
  xdouble *crl, *trl, **cr = NULL, **tr = NULL, **ck, **tk;
  xdouble **yr = NULL, *yrl = NULL, *arr, *vc = NULL;
  xdouble *ri, *ki, **r2p, **invr2p, **k2p, **invk2p, *rDm1, *kDm1;
  int i, dm, l, l0 = 1;
  long *coef;
  FFTWPFX(plan) plans[2] = {NULL, NULL};
  clock_t t0 = clock(), t1;

  rmax = adjustrmax(rmax, npt, &dr, &dm, ffttype);
  dk = PI/dr/npt;

  facr2k = pow_si(PI*2, K) *  dr;
  fack2r = pow_si(PI*2, -K-1) * dk;

  /* B2 = (PI*2)^K/(2 K + 1)!! */
  B2 = 1;
  for (i = 1; i <= K; i++)
    B2 *= PI*2/(2*i + 1);
  surfr = B2 * 2 * dim;
  surfk = surfr / pow_si(PI*2, dim);
  if ( gaussf ) B2 = SQRT(pow_si(PI, dim))/2;
  /* TODO: compute B2 for inverse potential */
  if ( invexp > 0 ) B2 = 0;
  printf("dr %f, dm %d, rmax %f, ffttype %d, HNC %d, B2 %g\n",
      (double) dr, dm, (double) rmax, ffttype, dohnc, (double) B2);

  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE2DARR(r2p, K, npt) /* r^{K - l} for l = 0, ..., K */
  MAKE2DARR(k2p, K, npt) /* k^{K - l} for l = 0, ..., K */
  MAKE2DARR(invr2p, K, npt) /* r^{-K-l} for l = 0, ..., K */
  MAKE2DARR(invk2p, K, npt) /* k^{-K-l} for l = 0, ..., K */
  MAKE1DARR(rDm1, npt); /* r^(D - 1) dr */
  MAKE1DARR(kDm1, npt); /* k^(D - 1) dk */
  {
    xdouble rl, invrl, kl, invkl;

    for ( i = 0; i < npt; i++ ) {
      ri[i]  = dr * (i*2 + (ffttype ? 1 : 0))/2;
      ki[i]  = dk * (i*2 + (ffttype ? 1 : 0))/2;
      rl = 1;
      kl = 1;
      for ( l = 1; l <= K; l++ ) {
        r2p[K - l][i] = (rl *= ri[i]);
        k2p[K - l][i] = (kl *= ki[i]);
      }

      if ( ffttype == 0 && i == 0 ) continue;

      invrl = 1/rl;
      invkl = 1/kl;
      for ( l = 0; l < K; l++ ) {
        invr2p[l][i] = invrl;
        invk2p[l][i] = invkl;
        invrl /= ri[i];
        invkl /= ki[i];
      }

      rDm1[i] = surfr * pow_si(ri[i], dim - 1) * dr;
      kDm1[i] = surfk * pow_si(ki[i], dim - 1) * dk;
    }
  }

  /* compute the coefficients of the spherical Bessel function */
  MAKE1DARR(coef, K);
  getjn(coef, K - 1);

  /* auxiliary array for FFTW
   * needs npt + 1 elements for FFTW_REDFT00 */
  MAKE1DARR(arr, npt + 1);

#ifndef NOFFTW
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
  MAKE2DARR(ck, nmax - 1, npt);
  MAKE2DARR(tk, nmax - 1, npt);

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
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( dohnc || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) {
      swr[i] = 1;
      if ( hncalpha >= 0 )
        swr[i] = 1 - EXP( -hncalpha * ri[i] );
    }
    for ( i = 0; i < npt; i++ ) yr[0][i] = hncamp / swr[i];
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  t1 = clock();
  if ( snapshot )
    l0 = snapshot_open(dim, nmax, rmax, dohnc, mkcorr, ring, singer,
        npt, ck, tk, cr, tr, crl, trl, yr);
  fnvir = savevirheadx(fnvir, systitle, dim, l0, nmax,
      dohnc, mkcorr, npt, rmax, t1 - t0,
      hncamp, hncq, hncalpha, shift);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* c_l(r) --> c_l(k) for the previous l */
    sphr(npt, crl, ck[l-1], facr2k, plans, arr, coef, r2p, invk2p, ffttype);

    if ( ring ) {
      /* compute the ring sum based on ck */
      Bh = get_ksum(l, npt, ck, kDm1, &Br);
      Br = (dohnc ? -Br * (l+1) : -Br * 2) / l;
    }

    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, plans, arr, coef, k2p, invr2p, ffttype);

    if ( tr != NULL ) {
      COPY1DARR(tr[l], trl, npt);
    }

    if ( yr != NULL ) /* compute the cavity function y(r) */
      get_yr_hncx(l, nmax, npt, yr, trl, hncq, swr);

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yr[l][i] - hncamp * hncq * trl[i];
    }

    if ( dohnc ) {
      /* hypernetted-chain approximation:
       * c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = (fr[i] + 1) * (yr[l][i] + (1 - hncamp * hncq) * trl[i]) - trl[i];
    } else {
      /* Percus-Yevick approximation:
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = fr[i] * trl[i];
    }

    /* in the PY case, y(r) = 1 + t(r) */
    yrl = dohnc ? yr[l] : trl;
    if ( smoothpot ) {
      Bv = get_Bv_hnc(npt, yrl, trl, hncamp, hncq, rdfr, rDm1, dim);
    } else {
      Bv = get_Bv_hnc_hs(yrl, trl, hncamp, hncq, dm, B2);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) surfr r^(D-1) dr */
    Bc = -integr(npt, crl, rDm1) / (l + 2);

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( dohnc ) {
        Bm = get_Bm_singer(l, npt, cr, tr, rDm1);
        Bh = get_Bh_singer(l, npt, cr, tr, rDm1) - Bh*(l+1)/2;
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, rDm1);
        Bh = get_Bp_py(l, npt, cr, tr, rDm1) - Bh;
      }
    } else {
      Bm = Bh = 0;
    }

    if ( mkcorr ) {
      /* for the Gaussian model, Bv is exact for l <= 2 (n <= 4)
       * otherwise, Bv is exact for l = 1 (n = 3) */
      if (l > 1 + gaussf) Bv *= 1 + shift;
      Bm = get_corr1x(l, npt, dm, crl, fr, rdfr, rDm1,
                      dim, B2, vc, &Bc, &Bv, &fcorr);
    }

    savevir(fnvir, dim, l+2, Bc, Bv, Bm, Bh, Br, B2, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yr);
    if ( snapshot )
      snapshot_take(l, npt, ck[l-1], tk[l], crl, trl, nmax, yr);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(coef, K);
#ifndef NOFFTW
  FFTWPFX(destroy_plan)(plans[0]);
  FFTWPFX(destroy_plan)(plans[1]);
#endif
  FREE1DARR(arr, npt);
  FREE1DARR(ri, npt);
  FREE1DARR(ki, npt);
  FREE2DARR(r2p, K, npt); FREE2DARR(invr2p, K, npt);
  FREE2DARR(k2p, K, npt); FREE2DARR(invk2p, K, npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(swr, npt);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(vc, npt);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, rmax, ffttype, dohnc);
#ifndef NOFFTW
  FFTWPFX(cleanup)();
#endif
  return 0;
}
