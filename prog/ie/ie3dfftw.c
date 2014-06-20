/* Computing the virial coefficients of the 3D hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ie3dfftw.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ie3dfftw.c -lfftw3l
 * Or for the 128-bit precision
 *  gcc -DF128 ieodfftw.c -lfftw3q -lquadmath -lm
 * To disable FFTW
 *  gcc -DNOFFTW ie3dfftw.c -lm
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



int nmax = 10;
xdouble rmax = 0;
int numpt = 1024;
int ffttype = 1;
int doHNC = 0;
int ring = 0;
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
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax <= 0 ) rmax = nmax + 2;
  if ( mkcorr ) singer = ring = 0;
  if ( singer ) ring = 1;
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*fac/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, xdouble *in, xdouble *out, xdouble fac,
    FFTWPFX(plan) p, xdouble *arr, xdouble *ri, xdouble *ki, int ffttype)
{
  int i, im = ffttype ? 0 : 1;

  for ( i = im; i < npt; i++ ) /* form in(x) * x */
    arr[i] = in[i] * ri[i];
#ifdef NOFFTW
  if ( ffttype )
    sint11(arr, npt);
  else
    sint00(arr, npt);
#else
  FFTWPFX(execute)(p);
#endif
  for ( i = im; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / ki[i];
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype, int doHNC)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk;
  xdouble Bc, Bv, Bm = 0, Bh, Br, B2, B2p, fcorr = 0;
  xdouble *fr, *crl, *trl, **ck, **tk, **cr = NULL, **tr = NULL;
  xdouble **yr = NULL, *arr, *ri, *ki, *ri2, *ki2, *vc = NULL;
  int i, dm, l;
  FFTWPFX(plan) plan = NULL;
  clock_t t0 = clock(), t1;

  rmax = adjustrmax(rmax, npt, &dr, &dm, ffttype);
  dk = PI/dr/npt;

  facr2k = PI*2 * dr;
  fack2r = pow_si(PI*2, -2) * dk;

  B2 = PI*2/3;
  surfr = PI*4;
  surfk = surfr * pow_si(PI*2, -3);
  printf("dr %f, dm %d, rmax %f, ffttype %d, HNC %d, B2 %g\n",
      (double) dr, dm, (double) rmax, ffttype, doHNC, (double) B2);

  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE1DARR(ri2, npt);
  MAKE1DARR(ki2, npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i]  = dr * (i*2 + (ffttype ? 1 : 0))/2;
    ki[i]  = dk * (i*2 + (ffttype ? 1 : 0))/2;
    ri2[i] = surfr * ri[i] * ri[i] * dr;
    ki2[i] = surfk * ki[i] * ki[i] * dk;
  }

  /* auxiliary array for FFTW */
  MAKE1DARR(arr, npt + 1);

#ifndef NOFFTW
  if ( ffttype ) {
    plan = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);
  } else {
    /* we only use npt - 1 points in this case */
    plan = FFTWPFX(plan_r2r_1d)(npt - 1, arr + 1, arr + 1, FFTW_RODFT00, FFTW_ESTIMATE);
  }
#endif

  MAKE1DARR(fr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE2DARR(ck, nmax - 1, npt);
  MAKE2DARR(tk, nmax - 1, npt);

  /* compute f(r) and f(k) = c0(k) */
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    crl[i] = fr[i] = (i < dm) ? -1 : 0;

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt);
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( doHNC || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) yr[0][i] = 1;
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  t1 = clock();
  fnvir = savevirhead(fnvir, NULL, 3, 1, nmax,
      doHNC, mkcorr, npt, rmax, t1 - t0);

  B2p = B2;
  for ( l = 1; l < nmax - 1; l++ ) { /* c_l and t_l, B_{l+2} */
    /* c_l(r) --> c_l(k) for the previous l */
    sphr(npt, crl, ck[l-1], facr2k, plan, arr, ri, ki, ffttype);

    if ( ring ) {
      /* compute the ring sum based on ck */
      Bh = get_ksum(l, npt, ck, ki2, &Br);
      Br = (doHNC ? -Br * (l+1) : -Br * 2) / l;
    }

    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, plan, arr, ki, ri, ffttype);

    if ( tr != NULL ) {
      COPY1DARR(tr[l], trl, npt);
    }

    if ( yr != NULL ) /* compute the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yr[l][i] - trl[i];
    }

    if ( doHNC ) {
      /* hypernetted-chain approximation:
       * c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = (fr[i] + 1) * yr[l][i] - trl[i];
      Bv = contactv(yr[l], dm, B2);
    } else {
      /* Percus-Yevick approximation:
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = fr[i] * trl[i];
      Bv = contactv(trl, dm, B2);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) surfr r^2 dr */
    Bc = -integr(npt, crl, ri2) / (l + 2);

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( doHNC ) {
        Bm = get_Bm_singer(l, npt, cr, tr, ri2);
        Bh = get_Bh_singer(l, npt, cr, tr, ri2) - Bh*(l+1)/2;
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, ri2);
        Bh = get_Bp_py(l, npt, cr, tr, ri2) - Bh;
      }
    } else {
      Bm = Bh = 0;
    }

    if ( mkcorr ) {
      /* in the PY case, y(r) = 1 + t(r) */
      Bm = get_corr1_hs(l, npt, dm, doHNC ? yr[l] : trl,
          crl, fr, ri2, B2, vc, &Bc, &Bv, &fcorr);
    }

    B2p *= B2;
    savevir(fnvir, 3, l+2, Bc, Bv, Bm, Bh, Br, B2p, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yr);
  }
  savevirtail(fnvir, clock() - t1);

#ifndef NOFFTW
  FFTWPFX(destroy_plan)(plan);
#endif
  FREE1DARR(arr, npt);
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
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, rmax, ffttype, doHNC);
#ifndef NOFFTW
  FFTWPFX(cleanup)();
#endif
  return 0;
}
