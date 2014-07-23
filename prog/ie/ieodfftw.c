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

int ietype = 0;

xdouble hncamp = 0, hncq = 1; /* Marucho-Pettitt */
xdouble hncalpha = -1; /* Rogers-Young */
xdouble hcs = 1; /* Hutchinson-Conkie s */
xdouble rowphi = 0; /* Rowlinson's Phi */
xdouble invphi = 0; /* inverse Rowlinson's Phi */
xdouble verleta = 0, verletb = 0; /* Verlet modified */
xdouble bbpgs = 15./8; /* MS/BBPG s */

xdouble shift = 0, shiftinc = 0;
int shiftl0 = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  xdouble shiftn = 0;
  int usems = 0;

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &dim, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", "%lf", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "-T", "%d", &ietype, "type of closure, 0: PY, 1: HNC-like; set the respective parameters automatically sets the type");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted-chain (HNC) like approximation");
  argopt_add(ao, "-q", "%" XDBLSCNF "f", &hncq,  "q of the hypernetted-chain (Marucho-Pettitt) approximation, y(r) = a0 exp(q t(r))");
  argopt_add(ao, "-a", "%" XDBLSCNF "f", &hncamp, "a0 of the hypernetted-chain (Marucho-Pettitt) approximation, 0: to be set by 1/q");
  argopt_add(ao, "--rya", "%" XDBLSCNF "f", &hncalpha, "alpha, in the Roger-Young switch function [1 - exp(-alpha*r)]^m");
  argopt_add(ao, "--hcs", "%" XDBLSCNF "f", &hcs, "s, in the Hutchinson-Conkie closure, y(r) = (1 + s t(r))^(1/s)");
  argopt_add(ao, "--rphi", "%" XDBLSCNF "f", &rowphi, "phi of the Rowlinson approximation, t(r) = (1 - phi) (y(r) - 1) + phi ln y(r)");
  argopt_add(ao, "--iphi", "%" XDBLSCNF "f", &invphi, "phi of the inverse Rowlinson approximation, y(r) = exp(t(r)) phi + (1 - phi) (1 + t(r))");
  argopt_add(ao, "--va", "%" XDBLSCNF "f", &verleta, "A of the Verlet approximation, y(r) = exp[t(r) - A t(r)^2 /2 / (1 + B t(r) / 2) ]");
  argopt_add(ao, "--vb", "%" XDBLSCNF "f", &verletb, "B of the Verlet approximation, y(r) = exp[t(r) - A t(r)^2 /2 / (1 + B t(r) / 2) ]");
  argopt_add(ao, "--bbpgs", "%" XDBLSCNF "f", &bbpgs, "s of the BBPG approximation, y(r) = exp[ (1 + s t(r))^(1/s) - 1 ]");
  argopt_add(ao, "--ms", "%b", &usems, "Martynov-Sarkisov approximation, y(r) = exp[ (1 + 2 t(r))^(1/2) - 1 ]");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "-c", "%" XDBLSCNF "f", &shift, "shift of t(r) in computing Bv");
  argopt_add(ao, "-C", "%" XDBLSCNF "f", &shiftn, "shift of t(r) in computing Bv, in terms of n");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &shiftinc, "increment of the shift");
  argopt_add(ao, "-L", "%d", &shiftl0, "minimal order l for the shift");
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
  if ( argopt_isset(ao, verleta) || argopt_isset(ao, verletb) ) ietype = IETYPE_VERLET;
  if ( argopt_isset(ao, usems) ) { ietype = IETYPE_BBPG; bbpgs = 2; }
  if ( argopt_isset(ao, bbpgs) ) ietype = IETYPE_BBPG;

  if ( ietype > IETYPE_HNC ) dohnc = 0;

  if ( dohnc ) ietype = IETYPE_HNC;

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
  if ( gaussf ) strcpy(systitle, "GF");
  else if ( invexp ) sprintf(systitle, "INV%d", invexp);

  printf("D %d, rmax %f, npt %d, ietype %d, HNC %d, a0 %g, q %g, %s\n",
      dim, (double) rmax, numpt, ietype, dohnc,
      (double) hncamp, (double) hncq, systitle);

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
  xdouble **yr = NULL, *yrl, *arr, *vc = NULL, *yrcoef = NULL;
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
  printf("D %d, dr %f, dm %d, rmax %f, ffttype %d, ietype %d, B2 %.10e\n",
      dim, (double) dr, dm, (double) rmax, ffttype, ietype, (double) B2);

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
  MAKE1DARR(yrl, npt);
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
  }

  if ( ietype > IETYPE_HNC || singer ) {
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( ietype > IETYPE_HNC ) { /* coefficients */
    MAKE1DARR(yrcoef, nmax - 1);
    if ( ietype == IETYPE_HC ) {
      init_hccoef(yrcoef, nmax - 1, hcs);
    } else if ( ietype == IETYPE_BBPG ) {
      init_bbpgcoef(yrcoef, nmax - 1, bbpgs);
    } else if ( ietype == IETYPE_ROWLINSON ) {
      init_rowlinsoncoef(yrcoef, nmax - 1, rowphi);
    } else if ( ietype == IETYPE_INVROWLINSON ) {
      init_invrowlinsoncoef(yrcoef, nmax - 1, invphi);
    } else if ( ietype == IETYPE_VERLET ) {
      init_verletcoef(yrcoef, nmax - 1, verleta, verletb);
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

  if ( dohnc || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) yr[0][i] = hncamp / swr[i];
  }

  if ( mkcorr ) {
    MAKE1DARR(vc, npt);
  }

  t1 = clock();

  if ( snapshot )
    l0 = snapshot_open(dim, nmax, rmax, dohnc, mkcorr, ring, singer,
        npt, ck, tk, cr, tr, crl, trl, yr);

  fnvir = savevirheadx(fnvir, systitle, dim, l0, nmax,
      ietype, mkcorr, npt, rmax, t1 - t0,
      hncamp, hncq, hncalpha, shift, shiftinc, shiftl0);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* c_l(r) --> c_l(k) for the previous l */
    sphr(npt, crl, ck[l-1], facr2k, plans, arr, coef, r2p, invk2p, ffttype);

    if ( ring ) {
      /* compute the ring sum based on ck */
      Bh = get_ksum(l, npt, ck, kDm1, &Br);
      Br = (dohnc ? -Br * (l+1) : -Br * 2) / l;
    }

    /* compute t_l(k) from c_0(k), ... c_{l-1}(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, plans, arr, coef, k2p, invr2p, ffttype);

    if ( tr != NULL ) { /* save tr if needed */
      COPY1DARR(tr[l], trl, npt);
    }

    /* compute the cavity function y(r) */
    if ( dohnc ) { /* specialized HNC closure */
      get_yr_hncx(l, nmax, npt, yr, trl, hncq, swr);
      for ( i = 0; i < npt; i++ )
        yrl[i] = yr[l][i] + (1 - hncamp * hncq) * trl[i];
    } else if ( ietype == IETYPE_PY) { /* PY closure */
      for ( i = 0; i < npt; i++ )
        yrl[i] = trl[i];
    } else {
      get_yr_series(l, npt, tr, yrcoef, yrl);
    }

    for ( i = 0; i < npt; i++ ) /* c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      crl[i] = (fr[i] + 1) * (mkcorr ? trl[i] : yrl[i]) - trl[i];

    Bv = get_Bv(npt, mkcorr ? trl : yrl, smoothpot, rdfr, rDm1, dim, dm, B2);
    Bc = get_Bc(l, npt, crl, rDm1);

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
      /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yrl[i] - trl[i];

      /* apply the shift */
      if ( l >= shiftl0 ) Bv *= 1 + shift + shiftinc * l;

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
  FREE1DARR(arr,  npt);
  FREE1DARR(ri,   npt);
  FREE1DARR(ki,   npt);
  FREE2DARR(r2p, K, npt); FREE2DARR(invr2p, K, npt);
  FREE2DARR(k2p, K, npt); FREE2DARR(invk2p, K, npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr,   npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(swr,  npt);
  FREE1DARR(crl,  npt);
  FREE1DARR(trl,  npt);
  FREE1DARR(yrl,  npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(vc, npt);
  FREE1DARR(yrcoef, nmax - 1);
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
