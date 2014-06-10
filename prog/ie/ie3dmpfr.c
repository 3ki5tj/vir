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
int singer = 0;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;



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
  if ( mkcorr ) singer = 0;
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
  mpfr_t dr, dk, pi2, facr2k, fack2r, surfr, surfk, B2, B2p, tmp1, tmp2, tmp3;
  mpfr_t Bc, Bv, Bm, Bh, Br, Bc0, dBc, Bv0, dBv, fcorr;
  mpfr_t *fr, *crl, *trl, **cr = NULL, **tr = NULL, **ck, **tk;
  mpfr_t **yr = NULL, *vc = NULL, *arr;
  mpfr_t *ri, *ki, *ri2, *ki2;
  double rmax;
  int i, dm, l;
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
  INIT_(tmp3);

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
  printf("precision %d, rmax %g, dk %g, %d bins in the hard core\n",
      i, rmax, GET_D_(dk), dm);

  MAKE1DARR(arr, npt + 1);

  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE1DARR(ri2, npt);
  MAKE1DARR(ki2, npt);

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

  /* surfr = PI*4 */
  MUL_SI_(surfr, pi2, 2);
  /* surfk = surfr / (2*PI)^3 */
  POW_SI_(tmp1, pi2, -3);
  MUL_(surfk, surfr, tmp1);
  for ( i = 0; i < npt; i++ ) {
    /* ri2[i] = ri[i]^2 * surfr * dr */
    SQR_(ri2[i], ri[i]);
    MUL3_X_(ri2[i], surfr, dr);
    /* ki2[i] = ki[i]^2 * surfk * dk */
    SQR_(ki2[i], ki[i]);
    MUL3_X_(ki2[i], surfk, dk);
  }

  /* B2 = PI*2/3; */
  DIV_SI_(B2, pi2, 3);

  MAKE1DARR(fr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE2DARR(ck, nmax - 1, npt)
  MAKE2DARR(tk, nmax - 1, npt)

  /* construct f(r) and f(k) */
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

  t1 = clock();
  fnvir = savevirhead(fnvir, NULL, 3, nmax, doHNC, mkcorr, npt, rmax, t1 - t0);

  SET_(B2p, B2);
  for ( l = 1; l < nmax - 1; l++ ) {
    if ( !mkcorr ) {
      /* compute the ring sum based on ck */
      get_ksum(Bh, l, npt, ck, ki2, Br);
      /* Br = (doHNC ? -Br * (l+1) : -Br * 2) / l; */
      MUL_SI_X_(Br, (doHNC ? -(l+1) : -2));
      DIV_SI_X_(Br, l);
    }

    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, arr, ki, ri, ffttype);

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

    if ( doHNC ) {
      /* hypernetted chain approximation: c(r) = (1 + f(r)) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ ) {
        ADD_SI_(tmp1, fr[i], 1);
        FMS_(crl[i], tmp1, yr[l][i], trl[i]);
      }
      /* Bv = B2*(yr[l][dm] + yr[l][dm-1])/2; */
      contactv(Bv, yr[l], dm, B2);
    } else {
      /* Percus-Yevick approximation: c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        MUL_(crl[i], fr[i], trl[i]);
      /* Bv = B2 * trl(1+) */
      contactv(Bv, trl, dm, B2);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    integr(Bc, npt, crl, ri2);
    DIV_SI_X_(Bc, -(l+2));

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( doHNC ) {
        get_Bm_singer(Bm, l, npt, cr, tr, ri2);
        get_Bh_singer(tmp1, l, npt, cr, tr, ri2);
        MUL_SI_(tmp2, Bh, l+1);
        DIV_SI_X_(tmp2, 2);
        SUB_(Bh, tmp1, tmp2);
        /* Bh = tmp1 - Bh*(l+1)/2; */
      } else {
        get_Bx_py(Bm, l, npt, cr, tr, ri2);
        get_Bp_py(tmp1, l, npt, cr, tr, ri2);
        SUB_(Bh, tmp1, Bh);
      }
    } else {
      SET_SI_(Bm, 0);
      SET_SI_(Bh, 0);
    }

    if ( mkcorr ) {
      get_corr1_hs(Bm, l, npt, dm, doHNC ? yr[l] : trl,
          crl, fr, ri2, B2, vc, Bc, Bv, fcorr);
    }

    MUL_X_(B2p, B2);
    savevir(fnvir, 3, l+2, Bc, Bv, Bm, Bh, Br, B2p, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yr);

    /* c_l(r) --> c_l(k) */
    sphr(npt, crl, ck[l], facr2k, arr, ri, ki, ffttype);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(arr, npt + 1);
  FREE1DARR(ri, npt);
  FREE1DARR(ki, npt);
  FREE1DARR(ri2, npt);
  FREE1DARR(ki2, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  if ( cr != NULL ) FREE2DARR(cr, nmax - 1, npt);
  if ( tr != NULL ) FREE2DARR(tr, nmax - 1, npt);
  if ( yr != NULL ) FREE2DARR(yr, nmax - 1, npt);
  if ( vc != NULL ) FREE1DARR(vc, npt);

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
  CLEAR_(tmp3);
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

  intgeq(nmax, numpt, rmax, ffttype, doHNC);

  MPFFT_ARR1D_FREE();
  mpfr_free_cache();
  ssdelall();
  return 0;
}
