/* Computing the virial coefficients of the Lennard-Jones fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ljie3d.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ljie3d.c -lfftw3l
 * To disable FFTW
 *  gcc -DNOFFTW ljie3dfftw.c -lm
 * */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#ifdef LDBL
typedef long double xdouble;
#define FFTWPFX(f) fftwl_##f
#define DBLSCNF "L"
#define DBLPRNF "L"
#else
typedef double xdouble;
#define FFTWPFX(f) fftw_##f
#define DBLSCNF "l"
#define DBLPRNF ""
#endif



#define PI (xdouble) 3.1415926535897932384626433832795L



#ifdef NOFFTW
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#else
#include <fftw3.h>
#endif

#include "ieutil.h"



int nmax = 10;
xdouble rmax = (xdouble) 10.24L;
int numpt = 1024;
xdouble beta = 1;
int ffttype = 1;
int doHNC = 0;
int singer = 0;
int mkcorr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  xdouble T = 1;
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-T", "%" DBLSCNF "f", &T, "temperature");
  argopt_add(ao, "-R", "%" DBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  printf("rmax = %f, T %lf, HNC %d\n", (double) rmax, (double) T, doHNC);
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



/* return the potential phi(r), and -r*phi'(r)*/
static xdouble pot(xdouble r, xdouble *ndphir)
{
  xdouble invr6 = 1/(r*r);
  invr6 = invr6*invr6*invr6;
  *ndphir = invr6*(48*invr6 - 24);
  return 4*invr6*(invr6 - 1);
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype, int doHNC)
{
  xdouble dr, dk, facr2k, fack2r, surf, B2, B2tail;
  xdouble *fr, *dfr, *crl, *trl, **ck, **tk, *Bc, *Bv;
  xdouble **yr = NULL, *arr, *ri, *ki, *ri2, *vc = NULL;
  xdouble **cr, **tr;
  int i, dm, l;
  FFTWPFX(plan) plan = NULL;

  dr = rmax / npt;
  dm = (int) (1/dr + .5); /* number of bins in the hard core */
  /* fix dr such that r = 1 lies at the middle of bins dm-1 and dm */
  if ( ffttype ) { /* r(i) = dr * (i + .5); */
    dr = (xdouble) 1/dm;
    rmax = dr * npt;
  } else { /* r(i) = dr * i; */
    dm += 1;
    dr = (xdouble) 1/(dm - .5);
    rmax = dr * (npt - .5);
  }
  dk = PI/dr/npt;

  MAKE1DARR(arr, npt);
  MAKE1DARR(ri, npt);
  MAKE1DARR(ki, npt);
  MAKE1DARR(ri2, npt);

#ifndef NOFFTW
  if ( ffttype ) {
    plan = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);
  } else {
    /* we only use npt - 1 points in this case */
    plan = FFTWPFX(plan_r2r_1d)(npt - 1, arr + 1, arr + 1, FFTW_RODFT00, FFTW_ESTIMATE);
  }
#endif
  for ( i = 0; i < npt; i++ ) {
    ri[i] = dr * (i*2 + (ffttype ? 1 : 0))/2;
    ki[i] = dk * (i*2 + (ffttype ? 1 : 0))/2;
  }

  facr2k = PI*2 * dr;
  fack2r = pow(PI*2, -2) * dk;

  MAKE1DARR(fr, npt);
  MAKE1DARR(dfr, npt);
  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  MAKE1DARR(trl, npt);
  MAKE1DARR(crl, npt);

  for ( i = 0; i < npt; i++ ) { /* compute f(r) = exp(-beta u(r)) - 1 */
    xdouble bphi = beta * pot(ri[i], &dfr[i]);
    fr[i] = exp(-bphi) - 1;
    dfr[i] *= beta * (1 + fr[i]);
  }
  sphr(npt, fr, ck[0], facr2k, plan, arr, ri, ki, ffttype); /* f(r) --> f(k) */

  surf = PI*4;
  for ( i = 0; i < npt; i++ )
    ri2[i] = surf * ri[i] * ri[i] * dr;
  B2 = -integr(npt, fr, ri2)/2;
  B2tail = -PI*8*pow(rmax, -3)/3;
  B2 += B2tail;
  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  Bc[2] = Bv[2] = B2;

  printf("beta %g, dr %f, dm %d, rmax %f, ffttype %d, HNC %d, B2 %g %g\n",
      (double) beta, (double) dr, dm, (double) rmax, ffttype, doHNC,
      (double) B2, (double) B2tail);

  if ( singer || 1 ) {
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

  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(npt, tk[l], trl, fack2r, plan, arr, ki, ri, ffttype);

    if ( yr != NULL ) { /* comput the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yr[l][i] - trl[i];
    }

    if ( doHNC ) {
      /* hypernetted approximation: c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++)
        crl[i] = (fr[i] + 1) * yr[l][i] - trl[i];
      Bv[l+2] = integr2(npt, dfr, yr[l], ri2) / 6;
    } else {
      /* Percus-Yevick approximation: c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = fr[i] * trl[i];
      Bv[l+2] = integr2(npt, dfr, trl, ri2) / 6;
    }

    if ( mkcorr ) {
      get_corr1_hnc(l, npt, doHNC ? yr[l] : trl, crl, fr, dfr, ri2, 3, vc);
      Bv[l+2] += integr2(npt, dfr, vc, ri2) / 6;
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    Bc[l+2] = -integr(npt, crl, ri2) / (l + 2);
    printf("Bc%-4d= %20.12" DBLPRNF "e, Bv%-4d= %20.12" DBLPRNF "e\n",
           l+2, Bc[l+2], l+2, Bv[l+2]);
    /* c(r) --> c(k) */
    sphr(npt, crl, ck[l], facr2k, plan, arr, ri, ki, ffttype);
  }

#ifndef NOFFTW
  FFTWPFX(destroy_plan)(plan);
#endif
  FREE1DARR(arr, npt);
  FREE1DARR(ri, npt);
  FREE1DARR(ki, npt);
  FREE1DARR(ri2, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(dfr, npt);
  FREE1DARR(trl, npt);
  FREE1DARR(crl, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE1DARR(Bc, nmax + 1);
  FREE1DARR(Bv, nmax + 1);
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

