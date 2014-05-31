/* Computing the virial coefficients of a odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ieodmp.c -lmpfr -lgmp
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpfft.h"
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
int nmax = 10;
char *rmax = NULL;
int numpt = 32768;
int ffttype = 1;
int doHNC = 0;
int mkcorr = 0;
int verbose = 0;



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
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax == NULL ) rmax = "32.768";
  if ( dim < 3 || dim % 2 == 0 ) argopt_help(ao);
  K = (dim - 1)/2;
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* get the coefficient c_l of spherical Bessel function jn(x)
 *  jn(x) = Sum_{l = 0 to n} c_l [sin(x) or cos(x)] / x^{n + l + 1},
 * where the l'th term is sin(x) if n + l is even, or cos(x) otherwise */
static void getjn(int *c, int n)
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
    printf(" %+d*%s/x^%d", c[i], fs[(i+n)%2], n + i + 1);
  printf("\n");
}



/* compute the D-dimensional Fourier transform of a spherical transform
 *    out(k) = 2 fac0 Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' is the intermediate array
 * */
static void sphr(int npt, mpfr_t *in, mpfr_t *out, mpfr_t fac,
    mpfr_t *arr, int *coef, mpfr_t **r2p, mpfr_t **k2q, int ffttype)
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



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, const char *srmax, int ffttype, int doHNC)
{
  mpfr_t dr, dk, pi2, facr2k, fack2r, surf, B2, tmp1, tmp2;
  mpfr_t Bc0, dBc, Bv0, dBv, eps;
  mpfr_t *fr, *crl, *trl, **ck, **tk, *Bc, *Bv;
  mpfr_t **yr = NULL, *arr, *vc;
  mpfr_t **r2p, **invr2p, **k2p, **invk2p, *rDm1;
  double rmax;
  int i, dm, l, *coef;

  INIT_(dr);
  INIT_(dk);
  INIT_(pi2);
  INIT_(facr2k);
  INIT_(fack2r);
  INIT_(surf);
  INIT_(B2);
  INIT_(tmp1);
  INIT_(tmp2);

  INIT_(Bc0);
  INIT_(dBc);
  INIT_(Bv0);
  INIT_(dBv);
  INIT_(eps);

  /* dr = rmax/npt */
  SET_STR_(dr, srmax);
  DIV_SI_(dr, dr, npt);
  /* dm: number of bins in the hard core */
  dm = (int) (1/GET_D_(dr) + .5);
  if ( ffttype ) {
    /* dr = 1./dm */
    SET_SI_(dr, 1);
    DIV_SI_X_(dr, dm);
    rmax = GET_D_(dr) * npt;
  } else {
    dm += 1;
    /* dr = 2./(dm*2 - 1) */
    SET_SI_(dr, 2);
    DIV_SI_X_(dr, dm*2 - 1);
    rmax = GET_D_(dr) * (npt - .5);
  }
  /* dk = PI/npt/dr */
  CONST_PI_(dk);
  DIV_SI_X_(dk, npt);
  DIV_X_(dk, dr);

  /* print out the basic information */
  i = (int) mpfr_get_default_prec();
  printf("precision %d, rmax %g, dk %g, %d bins in the hard core\n",
      i, rmax, GET_D_(dk), dm);

  /* compute the coefficients of the spherical Bessel function */
  xnew(coef, K);
  getjn(coef, K - 1);

  /* auxiliary array, needs npt + 1 elements for cosine type-0 */
  MAKE1DARR(arr, npt + 1);

  MAKE2DARR(r2p, K, npt) /* r^{K - l} for l = 0, ..., K */
  MAKE2DARR(k2p, K, npt) /* k^{K - l} for l = 0, ..., K */
  MAKE2DARR(invr2p, K, npt) /* r^{-K-l} for l = 0, ..., K */
  MAKE2DARR(invk2p, K, npt) /* k^{-K-l} for l = 0, ..., K */
  MAKE1DARR(rDm1, npt); /* r^(D - 1) dr */
  {
    mpfr_t ri, rl, invrl, ki, kl, invkl;

    INIT_(ri);
    INIT_(rl);
    INIT_(invrl);
    INIT_(ki);
    INIT_(kl);
    INIT_(invkl);
    for ( i = 0; i < npt; i++ ) {
      /* ri = dr * (i*2 + 1)/2; */
      MUL_SI_(ri, dr, i*2 + (ffttype ? 1 : 0));
      DIV_SI_X_(ri, 2);
      /* rl = 1 */
      SET_SI_(rl, 1);

      /* ki = dk * (i*2 + 1)/2; */
      MUL_SI_(ki, dk, i*2 + (ffttype ? 1 : 0));
      DIV_SI_X_(ki, 2);
      /* kl = 1 */
      SET_SI_(kl, 1);
      for ( l = 1; l <= K; l++ ) {
        MUL_X_(rl, ri); /* rl *= ri */
        SET_(r2p[K-l][i], rl); /* r2p[K-l][i] = rl; */
        MUL_X_(kl, ki); /* kl *= ki */
        SET_(k2p[K-l][i], kl); /* k2p[K-l][i] = kl; */
      }

      if ( ffttype == 0 && i == 0 ) continue;

      SI_DIV_(invrl, 1, rl); /* invrl = 1/rl; */
      SI_DIV_(invkl, 1, kl); /* invkl = 1/kl; */
      for ( l = 0; l < K; l++ ) {
        SET_(invr2p[l][i], invrl); /* invr2p[l][i] = invrl; */
        SET_(invk2p[l][i], invkl); /* invk2p[l][i] = invkl; */
        DIV_X_(invrl, ri); /* invrl /= ri; */
        DIV_X_(invkl, ki); /* invkl /= ki; */
      }
    }
    CLEAR_(ri);
    CLEAR_(rl);
    CLEAR_(invrl);
    CLEAR_(ki);
    CLEAR_(kl);
    CLEAR_(invkl);
  }

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
  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  /* Bc[2] = Bv[2] = B2; */
  SET_(Bc[2], B2);
  SET_(Bv[2], B2);

  MUL_SI_(surf, B2, dim*2);
  for ( i = 0; i < npt; i++ ) {
    /* tmp1 = surf * ri^(dim - 1) * dr */
    POW_SI_(tmp1, r2p[K-1][i], dim - 1);
    MUL3_(rDm1[i], tmp1, surf, dr);
  }

  MAKE1DARR(fr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)

  /* construct f(r) and f(k) */
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    SET_SI_(fr[i], (i < dm) ? -1 : 0);
  /* f(r) --> f(k) */
  sphr(npt, fr, ck[0], facr2k, arr, coef, r2p, invk2p, ffttype);
  //for (i = 0; i<npt;i++){printf("%d %g %g\n", i, GET_D_(k2p[K-1][i]), GET_D_(ck[0][i]));}// exit(1);

  if ( doHNC || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ )
      SET_SI_(yr[0][i], 1);
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, arr, coef, k2p, invr2p, ffttype);

    if ( yr != NULL ) { /* compute the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        SUB_(vc[i], yr[l][i], trl[i]);
    }

    if ( doHNC ) {
      /* hypernetted chain approximation:
       * c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ ) {
        ADD_SI_(tmp1, fr[i], 1);
        FMS_(crl[i], tmp1, yr[l][i], trl[i]);
      }
      /* Bv[l+2] = B2*(yr[l][dm] + yr[l][dm-1])/2; */
      contactv(Bv[l+2], yr[l], dm, B2);
    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        MUL_(crl[i], fr[i], trl[i]);
      /* Bv[l+2] = B2*(trl[dm] + trl[dm-1])/2; */
      contactv(Bv[l+2], trl, dm, B2);
    }

    if ( mkcorr ) {
      get_corr1_hnc_hs(l, npt, dm, doHNC ? yr[l] : trl, crl, fr, rDm1, B2, vc);
      contactv(tmp1, vc, dm, B2);
      mpfr_printf("Bv %Rg %Rg\n", Bv[l+2], tmp1);
      ADD_X_(Bv[l+2], tmp1);
      mpfr_printf("Bv %Rg %Rg\n", Bv[l+2], tmp1);
    }

    /* Bc_{l+2} = -[1/(l+2)] Int c_l(r) S_D r^(D-1) dr */
    integr(Bc[l+2], npt, crl, rDm1);
    DIV_SI_X_(Bc[l+2], -(l+2));

    /* tmp1 = Bc/B2^(l+1), tmp2 = Bv/B2^(l+1) */
    POW_SI_(tmp2, B2, l+1);
    DIV_(tmp1, Bc[l+2], tmp2);
    DIV_(tmp2, Bv[l+2], tmp2);
    mpfr_printf("Bc(%3d) = %16.9Re (%16.9Re), "
                "Bv(%3d) = %16.9Re (%16.9Re)\n",
                l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c_l(r) --> c_l(k) */
    sphr(npt, crl, ck[l], facr2k, arr, coef, r2p, invk2p, ffttype);
  }

  FREE1DARR(arr, npt + 1);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE1DARR(fr, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE1DARR(Bc, nmax + 1);
  FREE1DARR(Bv, nmax + 1);
  free(coef);
  FREE1DARR(rDm1, npt);
  FREE2DARR(r2p, K, npt); FREE2DARR(invr2p, K, npt);
  FREE2DARR(k2p, K, npt); FREE2DARR(invk2p, K, npt);
  if ( doHNC || mkcorr ) {
    FREE2DARR(yr, nmax - 1, npt);
    if ( mkcorr ) {
      FREE1DARR(vc, npt);
    }
  }

  CLEAR_(dr);
  CLEAR_(dk);
  CLEAR_(pi2);
  CLEAR_(facr2k);
  CLEAR_(fack2r);
  CLEAR_(surf);
  CLEAR_(B2);
  CLEAR_(tmp1);
  CLEAR_(tmp2);
  CLEAR_(Bc0);
  CLEAR_(dBc);
  CLEAR_(Bv0);
  CLEAR_(dBv);
  CLEAR_(eps);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mpfr_set_default_prec(prec);

  intgeq(nmax, numpt, rmax, ffttype, doHNC);

  MPFFT_ARR1D_FREE();
  mpfr_free_cache();
  ssdelall();
  return 0;
}