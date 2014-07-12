/* Computing the virial coefficients of a odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ieodmp.c -lmpfr -lgmp
 * */
#include <stdio.h>
#include <math.h>
#include "fftmpfr.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"

#include "ieutilmpfr.h"



#ifdef D
int dim = D;
#else
int dim = 3;
#endif

int K;
int prec = 256;
int nmax = 12;
char *rmax = NULL;
int numpt = 32768;
int ffttype = 1;
int dohnc = 0;
int ring = 0;
int singer = 0;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the PY/HNC closure for odd-dimensional hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", NULL, &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-p", "%d", &prec, "float-point precision in bits");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted chain approximation");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax == NULL ) {
    static char srmax[32];
    sprintf(rmax = srmax, "%d", nmax + 2);
  }
  if ( dim < 3 || dim % 2 == 0 ) argopt_help(ao);
  K = (dim - 1)/2;
  if ( mkcorr ) singer = ring = 0;
  if ( singer ) ring = 1;
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



/* compute the D-dimensional Fourier transform of a spherical transform
 *    out(k) = 2 fac0 Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' is the intermediate array
 * */
static void sphr(int npt, mpfr_t *in, mpfr_t *out, mpfr_t fac,
    mpfr_t *arr, long *coef, mpfr_t **r2p, mpfr_t **k2q, int ffttype)
{
  int i, l, iscos;
  mpfr_t x;

  INIT_(x);

  /* clear the output */
  for ( i = 0; i < npt; i++ )
    SET_SI_(out[i], 0); /* out[i] = 0; */

  /* several rounds of transforms */
  for ( l = 0; l < K; l++ ) {
    /* decide if we want to do a sine or cosine transform */
    iscos = (K + l + 1) % 2;
    for ( i = 0; i < npt; i++ ) {
      /*
       * arr[i] = in[i] * r^(K - l) * coef[l];
       * or
       * arr[i] = in[i] * r2p[l][i] * coef[l];
       * */
      MUL_SI_(x, r2p[l][i], coef[l]);
      MUL_(arr[i], in[i], x);
    }

    /* do the sine or cosine transform */
    if ( iscos ) {
      if ( ffttype ) {
        mpcost11(arr, npt, NULL);
      } else {
        mpcost00(arr, npt, NULL);
        SET_SI_(arr[npt], 0);
      }
    } else {
      if ( ffttype )
        mpsint11(arr, npt, NULL);
      else
        mpsint00(arr, npt, NULL);
    }

    for ( i = 0; i < npt; i++ ) {
      /* out[i] += arr[i] * fac / k^(K + l);
       * or
       * out[i] += arr[i] * fac * k2q[l][i];
       * */
      MUL3_(x, arr[i], fac, k2q[l][i]);
      ADD_X_(out[i], x);
    }
  }

  CLEAR_(x);
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, const char *srmax, int ffttype, int dohnc)
{
  mpfr_t dr, dk, pi2, facr2k, fack2r, surfr, surfk, B2, B2p, tmp1, tmp2;
  mpfr_t Bc, Bv, Bm, Bh, Br, Bc0, dBc, Bv0, dBv, fcorr;
  mpfr_t *fr, *crl, *trl, **cr = NULL, **tr = NULL, **ck, **tk;
  mpfr_t **yr = NULL, *vc = NULL, *arr;
  mpfr_t *ri, *ki, **r2p, **invr2p, **k2p, **invk2p, *rDm1, *kDm1;
  double rmax;
  int i, dm, l;
  long *coef;
  clock_t t0 = clock(), t1;

  INIT_(dr);
  INIT_(dk);
  INIT_(pi2);
  INIT_(facr2k);
  INIT_(fack2r);
  INIT_(surfr);
  INIT_(surfk);
  INIT_(B2);
  INIT_(B2p);
  INIT_(tmp1);
  INIT_(tmp2);

  INIT_(Bc);
  INIT_(Bv);
  INIT_(Bm);
  INIT_(Bh);
  INIT_(Br);
  INIT_(Bc0);
  INIT_(dBc);
  INIT_(Bv0);
  INIT_(dBv);
  INIT_(fcorr);

  SET_SI_(Bm, 0);
  SET_SI_(Bh, 0);
  SET_SI_(Br, 0);

  /* dm: number of bins in the hard core */
  rmax = adjustrmax(srmax, npt, dr, &dm, ffttype);
  /* dk = PI/npt/dr */
  CONST_PI_(dk);
  DIV_SI_X_(dk, npt);
  DIV_X_(dk, dr);

  /* print out the basic information */
  i = (int) mpfr_get_default_prec();
  printf("D %d, precision %d, rmax %g, dk %g, %d bins in the hard core\n",
      dim, i, rmax, GET_D_(dk), dm);

  /* pi2 = PI*2 */
  CONST_PI_(pi2);
  MUL_SI_X_(pi2, 2);

  /* facr2k = (PI*2)^K */
  POW_SI_(facr2k, pi2, K);
  MUL_X_(facr2k, dr); /* facr2k *= dr; */

  /* fack2r = (PI*2)^(-K-1) */
  POW_SI_(fack2r, pi2, -K-1);
  MUL_X_(fack2r, dk); /* fack2r *= dk */

  /* B2 = (PI*2)^K/(2 K + 1)!! */
  SET_SI_(B2, 1); /* B2 = 1 */
  for (i = 1; i <= K; i++) {
    /* B2 *= PI*2/(i*2 + 1); */
    MUL_X_(B2, pi2);
    DIV_SI_X_(B2, i*2 + 1);
  }

  MUL_SI_(surfr, B2, dim*2);
  POW_SI_(tmp1, pi2, dim);
  DIV_(surfk, surfr, tmp1);

  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE2DARR(r2p, K, npt) /* r^{K - l} for l = 0, ..., K */
  MAKE2DARR(k2p, K, npt) /* k^{K - l} for l = 0, ..., K */
  MAKE2DARR(invr2p, K, npt) /* r^{-K-l} for l = 0, ..., K */
  MAKE2DARR(invk2p, K, npt) /* k^{-K-l} for l = 0, ..., K */
  MAKE1DARR(rDm1, npt); /* r^(D - 1) dr */
  MAKE1DARR(kDm1, npt); /* k^(D - 1) dk */
  {
    mpfr_t rl, invrl, kl, invkl;

    INIT_(rl);
    INIT_(invrl);
    INIT_(kl);
    INIT_(invkl);
    for ( i = 0; i < npt; i++ ) {
      /* ri = dr * (i*2 + 1)/2; */
      MUL_SI_(ri[i], dr, i*2 + (ffttype ? 1 : 0));
      DIV_SI_X_(ri[i], 2);
      /* rl = 1 */
      SET_SI_(rl, 1);

      /* ki = dk * (i*2 + 1)/2; */
      MUL_SI_(ki[i], dk, i*2 + (ffttype ? 1 : 0));
      DIV_SI_X_(ki[i], 2);
      /* kl = 1 */
      SET_SI_(kl, 1);
      for ( l = 1; l <= K; l++ ) {
        MUL_X_(rl, ri[i]); /* rl *= ri */
        SET_(r2p[K-l][i], rl); /* r2p[K-l][i] = rl; */
        MUL_X_(kl, ki[i]); /* kl *= ki */
        SET_(k2p[K-l][i], kl); /* k2p[K-l][i] = kl; */
      }

      if ( ffttype == 0 && i == 0 ) continue;

      SI_DIV_(invrl, 1, rl); /* invrl = 1/rl; */
      SI_DIV_(invkl, 1, kl); /* invkl = 1/kl; */
      for ( l = 0; l < K; l++ ) {
        SET_(invr2p[l][i], invrl); /* invr2p[l][i] = invrl; */
        SET_(invk2p[l][i], invkl); /* invk2p[l][i] = invkl; */
        DIV_X_(invrl, ri[i]); /* invrl /= ri; */
        DIV_X_(invkl, ki[i]); /* invkl /= ki; */
      }

      /* tmp1 = surf * ri^(dim - 1) * dr */
      POW_SI_(tmp1, r2p[K-1][i], dim - 1);
      MUL3_(rDm1[i], tmp1, surfr, dr);
      POW_SI_(tmp1, k2p[K-1][i], dim - 1);
      MUL3_(kDm1[i], tmp1, surfk, dk);
    }
    CLEAR_(rl);
    CLEAR_(invrl);
    CLEAR_(kl);
    CLEAR_(invkl);
  }

  /* compute the coefficients of the spherical Bessel function */
  xnew(coef, K);
  getjn(coef, K - 1);

  /* auxiliary array for FFT
   * needs npt + 1 elements for cosine type-0 */
  MAKE1DARR(arr, npt + 1);

  MAKE1DARR(fr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE2DARR(ck, nmax - 1, npt);
  MAKE2DARR(tk, nmax - 1, npt);

  /* construct f(r) and f(k) */
  for ( i = 0; i < npt; i++ ) { /* compute f(r) = exp(-beta u(r)) - 1 */
    SET_SI_(fr[i], (i < dm) ? -1 : 0);
    SET_(crl[i],  fr[i]);
  }

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt);
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( dohnc || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ )
      SET_SI_(yr[0][i], 1);
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  t1 = clock();
  fnvir = savevirhead(fnvir, NULL, dim, nmax, dohnc, mkcorr, npt, rmax, t1 - t0);

  SET_(B2p, B2);
  for ( l = 1; l < nmax - 1; l++ ) {
    /* c_l(r) --> c_l(k) for the previous l */
    sphr(npt, crl, ck[l-1], facr2k, arr, coef, r2p, invk2p, ffttype);

    if ( ring ) {
      /* compute the ring sum based on ck */
      get_ksum(Bh, l, npt, ck, kDm1, Br);
      /* Br = (dohnc ? -Br * (l+1) : -Br * 2) / l; */
      MUL_SI_X_(Br, (dohnc ? -(l+1) : -2));
      DIV_SI_X_(Br, l);
    }

    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, arr, coef, k2p, invr2p, ffttype);

    if ( tr != NULL ) {
      COPY1DARR(tr[l], trl, npt);
    }

    if ( yr != NULL ) { /* compute the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        SUB_(vc[i], yr[l][i], trl[i]);
    }

    if ( dohnc ) {
      /* hypernetted chain approximation:
       * c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ ) {
        ADD_SI_(tmp1, fr[i], 1);
        FMS_(crl[i], tmp1, yr[l][i], trl[i]);
      }
      /* Bv = B2*(yr[l][dm] + yr[l][dm-1])/2; */
      contactv(Bv, yr[l], dm, B2);
    } else {
      /* Percus-Yevick approximation:
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        MUL_(crl[i], fr[i], trl[i]);
      /* Bv = B2*(trl[dm] + trl[dm-1])/2; */
      contactv(Bv, trl, dm, B2);
    }

    /* Bc_{l+2} = -[1/(l+2)] Int c_l(r) S_D r^(D-1) dr */
    integr(Bc, npt, crl, rDm1);
    DIV_SI_X_(Bc, -(l+2));

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( dohnc ) {
        get_Bm_singer(Bm, l, npt, cr, tr, rDm1);
        get_Bh_singer(tmp1, l, npt, cr, tr, rDm1);
        MUL_SI_(tmp2, Bh, l+1);
        DIV_SI_X_(tmp2, 2);
        SUB_(Bh, tmp1, tmp2);
        /* Bh = tmp1 - Bh*(l+1)/2; */
      } else {
        get_Bx_py(Bm, l, npt, cr, tr, rDm1);
        get_Bp_py(tmp1, l, npt, cr, tr, rDm1);
        SUB_(Bh, tmp1, Bh);
      }
    } else {
      SET_SI_(Bm, 0);
      SET_SI_(Bh, 0);
    }

    if ( mkcorr ) {
      get_corr1_hs(Bm, l, npt, dm, dohnc ? yr[l] : trl,
          crl, fr, rDm1, B2, vc, Bc, Bv, fcorr);
    }
    MUL_X_(B2p, B2);
    savevir(fnvir, dim, l+2, Bc, Bv, Bm, Bh, Br, B2p, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yr);
  }
  savevirtail(fnvir, clock() - t1);

  free(coef);
  FREE1DARR(arr, npt + 1);
  FREE1DARR(ri, npt);
  FREE1DARR(ki, npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE2DARR(r2p, K, npt); FREE2DARR(invr2p, K, npt);
  FREE2DARR(k2p, K, npt); FREE2DARR(invk2p, K, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(vc, npt);

  CLEAR_(dr);
  CLEAR_(dk);
  CLEAR_(pi2);
  CLEAR_(facr2k);
  CLEAR_(fack2r);
  CLEAR_(surfr);
  CLEAR_(surfk);
  CLEAR_(B2);
  CLEAR_(B2p);
  CLEAR_(tmp1);
  CLEAR_(tmp2);
  CLEAR_(Bc);
  CLEAR_(Bv);
  CLEAR_(Bm);
  CLEAR_(Bh);
  CLEAR_(Br);
  CLEAR_(Bc0);
  CLEAR_(dBc);
  CLEAR_(Bv0);
  CLEAR_(dBv);
  CLEAR_(fcorr);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mpfr_set_default_prec(prec);

  intgeq(nmax, numpt, rmax, ffttype, dohnc);

  MPFFT_ARR1D_FREE();
  mpfr_free_cache();
  ssdelall();
  return 0;
}
