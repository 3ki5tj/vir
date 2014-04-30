/* Computing the virial coefficients of the 3D hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ie3dmp.c -lmpfr -lgmp
 * */
#include <stdio.h>
#include <math.h>
#include "mpfft.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#include "ieutilmpfr.h"



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

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", NULL, &rmax, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-p", "%d", &prec, "float-point precision in bits");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax == NULL ) rmax = "32.768";
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, mpfr_t *in, mpfr_t *out,
    mpfr_t fac, mpfr_t *arr, mpfr_t *ri, mpfr_t *ki, int ffttype)
{
  int i, im = ffttype ? 0 : 1;

  for ( i = im; i < npt; i++ ) { /* form in(x) * x */
    /* arr[i] = in[i] * ri[i]; */
    MUL_(arr[i], in[i], ri[i]);
  }

  /* do the sine transform */
  if ( ffttype )
    mpsint11(arr, npt, NULL);
  else
    mpsint00(arr, npt, NULL);

  for ( i = im; i < npt; i++ ) { /* form out(k) / k */
    /* out[i] = arr[i] * fac / ki[i] */
    MUL_(out[i], arr[i], fac);
    DIV_X_(out[i], ki[i]);
  }
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, const char *srmax, int ffttype, int doHNC)
{
  mpfr_t dr, dk, pi2, facr2k, fack2r, surf, B2, tmp1, tmp2, tmp3;
  mpfr_t Bc0, dBc, Bv0, dBv, eps;
  mpfr_t **ck, **tk, *fr, *trl, *crl, *Bc, *Bv, *Bm;
  mpfr_t **yr = NULL, *arr, *ri, *ki, *ri2, *vc;
  double rmax;
  int i, dm, l;

  INIT_(dr);
  INIT_(dk);
  INIT_(pi2);
  INIT_(facr2k);
  INIT_(fack2r);
  INIT_(surf);
  INIT_(B2);
  INIT_(tmp1);
  INIT_(tmp2);
  INIT_(tmp3);

  INIT_(Bc0);
  INIT_(dBc);
  INIT_(Bv0);
  INIT_(dBv);
  INIT_(eps);

  /* dr = rmax/npt */
  SET_STR_(dr, srmax);
  DIV_SI_X_(dr, npt);
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

  MAKE1DARR(arr, npt);
  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE1DARR(ri2, npt);

  for ( i = 0; i < npt; i++ ) {
    MUL_SI_(ri[i], dr, i*2 + (ffttype ? 1 : 0));
    DIV_SI_X_(ri[i], 2);
    MUL_SI_(ki[i], dk, i*2 + (ffttype ? 1 : 0));
    DIV_SI_X_(ki[i], 2);
  }

  CONST_PI_(pi2);
  MUL_SI_X_(pi2, 2);

  MUL_(facr2k, pi2, dr); /* facr2k = 2*pi * dr */

  POW_SI_(fack2r, pi2, -2); /* fack2r = 1/(2*pi)^2 */
  MUL_X_(fack2r, dk); /* fack2r *= dk */

  /* surf = PI*4 */
  MUL_SI_(surf, pi2, 2);
  for ( i = 0; i < npt; i++ ) {
    /* ri2[i] = ri[i]^2 * surf * dr */
    SQR_(ri2[i], ri[i]);
    MUL3_X_(ri2[i], surf, dr);
  }

  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  MAKE1DARR(Bm, nmax + 1);
  /* B2 = PI*2/3; */
  DIV_SI_(B2, pi2, 3);
  SET_(Bv[2], B2);
  SET_(Bc[2], B2);
  SET_(Bm[2], B2);

  /* construct f(r) and f(k) */
  MAKE1DARR(fr, npt);
  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  MAKE1DARR(trl, npt);
  MAKE1DARR(crl, npt);

  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    SET_SI_(fr[i], (i < dm) ? -1 : 0);
  sphr(npt, fr, ck[0], facr2k, arr, ri, ki, ffttype); /* f(r) --> f(k) */

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
    sphr(npt, tk[l], trl, fack2r, arr, ki, ri, ffttype);

    if ( yr != NULL ) { /* compute the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        SUB_(vc[i], yr[l][i], trl[i]);
    }

    if ( doHNC ) {
      /* hypernetted chain approximation: c(r) = (1 + f(r)) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ ) {
        ADD_SI_(tmp1, fr[i], 1);
        FMS_(crl[i], tmp1, yr[l][i], trl[i]);
      }
      /* Bv[l+2] = B2*(yr[l][dm] + yr[l][dm-1])/2; */
      contactv(Bv[l+2], yr[l], dm, B2);
    } else {
      /* Percus-Yevick approximation: c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        MUL_(crl[i], fr[i], trl[i]);
      /* Bv(l+2) = B2 * trl(1+) */
      contactv(Bv[l+2], trl, dm, B2);
    }

    if ( mkcorr ) {
      get_corr1_hnc_hs(l, npt, dm, doHNC ? yr[l] : trl, crl, fr, ri2, B2, vc);
      contactv(tmp1, vc, dm, B2);
      ADD_X_(Bv[l+2], tmp1);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    integr(Bc[l+2], npt, crl, ri2);
    DIV_SI_X_(Bc[l+2], -(l+2));

    /* tmp1 = Bc/B2^(l+1), tmp2 = Bv/B2^(l+1) */
    POW_SI_(tmp3, B2, l+1);
    DIV_(tmp1, Bc[l+2], tmp3);
    DIV_(tmp2, Bv[l+2], tmp3);
    mpfr_printf("Bc(%3d) = %16.9Re (%16.9Re), "
                "Bv(%3d) = %16.9Re (%16.9Re)\n",
                l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c(r) --> c(k) */
    sphr(npt, crl, ck[l], facr2k, arr, ri, ki, ffttype);
  }

  FREE1DARR(arr, npt);
  FREE1DARR(ri, npt);
  FREE1DARR(ki, npt);
  FREE1DARR(ri2, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(trl, npt);
  FREE1DARR(crl, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE1DARR(Bc, nmax + 1);
  FREE1DARR(Bv, nmax + 1);
  FREE1DARR(Bm, nmax + 1);
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
  CLEAR_(tmp3);
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
