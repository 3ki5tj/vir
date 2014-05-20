/* free energy related quantities from the HNC closure */
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#include "xdouble.h"



#define PI (xdouble) 3.1415926535897932384626433832795L



#ifdef NOFFTW
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#else
#include <fftw3.h>
#endif



int numpt = 8192;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 10;
xdouble beta;
xdouble rho = (xdouble) 0.8L;
int itermax = 10000;
xdouble tol = (xdouble) 1e-12L;
xdouble delta = (xdouble) 0.0005L;
xdouble sigvv = 1;
xdouble siguv = 2;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing free energy related quantities from the HNC closure";
  argopt_add(ao, "-T", "%" DBLSCNF "f", &T, "temperature");
  argopt_add(ao, "-r", "%" DBLSCNF "f", &rho, "rho");
  argopt_add(ao, "-R", "%" DBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-d", "%" DBLSCNF "f", &delta, "delta lambda");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  printf("rmax %f, rho %g, T %f\n",
      (double) rmax, (double) rho, (double) T);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*fac/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, xdouble *in, xdouble *out, xdouble fac,
    FFTWPFX(plan) p, xdouble *arr, xdouble *ri, xdouble *ki)
{
  int i;

  for ( i = 0; i < npt; i++ ) /* form in(x) * x */
    arr[i] = in[i] * ri[i];
#ifdef NOFFTW
  sint11(arr, npt);
#else
  FFTWPFX(execute)(p);
#endif
  for ( i = 0; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / ki[i];
}



/* return the potential phi(r), and -r*phi'(r)*/
static xdouble pot(xdouble r, xdouble sig, xdouble eps, xdouble *ndphir)
{
  xdouble invr6 = (sig*sig)/(r*r), u;
  invr6 = invr6*invr6*invr6;
  u = 4*invr6*(invr6 - 1);
  if (u > 1000) u = 1000;
  *ndphir = invr6*(48*invr6 - 24);
  if (*ndphir > 12000) *ndphir = 12000;
  *ndphir *= eps;
  return eps*u;
}



static void iter(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *fr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itermax)
{
  int i, iter;
  xdouble x, err, errmax;

  for ( i = 0; i < npt; i++ ) cr[i] = fr[i];
  for ( iter = 0; iter < itermax; iter++ ) {
    sphr(npt, cr, ck, facr2k, plan, arr, ri, ki);
    for ( i = 0; i < npt; i++ ) {
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    }
    sphr(npt, tk, tr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = (1 + fr[i]) * exp(tr[i]) - 1 - tr[i];
      //x = (1 + fr[i]) * (1 + tr[i] + tr[i]*tr[i]/2) - 1 - tr[i];
      if ((err = fabs(cr[i] - x)) > errmax) errmax = err;
      cr[i] = x;
    }
    if ( errmax < 1e-8 ) break;
  }
  //printf("iter %d errmax %g\n", iter, (double) errmax);
}



static void output(int npt, xdouble *ri,
    xdouble *cr, xdouble *tr, xdouble *bphi, char *fn)
{
  int i;
  FILE *fp;

  xfopen(fp, fn, "w", return);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%8.6f %14.8f %14.8f %g\n", (double) ri[i],
        (double) cr[i], (double) tr[i], (double) bphi[i]);
  }
  fclose(fp);
}



static xdouble getfe_hnc(int npt, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ri2,
    xdouble *ck, xdouble *tk, xdouble *ki2,
    xdouble *mu, xdouble *compr, xdouble *ddmu)
{
  int i;
  xdouble fe, x;

  fe = 0;
  *mu = 0;
  *compr = 0;
  *ddmu = 0;
  for ( i = 0; i < npt; i++ ) {
    /* 2 beta F = - rho^2 Int [c(r) - h(r)^2/2] dr
     *            + Int [log(1 - rho c(k)) + rho c(k)] dk/(2pi)^3
     *          = - rho^2 Int [c(r) - t(r)^2/2 - t(r) c(r)] dr
     *            + Int [log(1 - rho c(k)) + rho c(k) + rho^2 c^2(k)/2] dk/(2pi)^3
     *          = - Sum_G I(G) / s(G) */
    fe -= (cr[i] - tr[i] * (tr[i]*.5 + cr[i])) * ri2[i];
    /* mu = beta mu = -rho^2 Int [c(r) - t(r) h(r)/2] dr
     *    = - Sum_G I(G) n(G) / s(G) */
    *mu -= (cr[i] - tr[i] * (cr[i] + tr[i])*.5) * ri2[i];
    /* compr = d(beta mu) / drho = -Int c(r) dr */
    *compr -= cr[i] * ri2[i];
    *ddmu -= tr[i]*(tr[i] + cr[i]) * ri2[i];
  }
  fe *= rho * rho;
  *mu *= rho;
  *ddmu /= rho;
  /* the k-space part of the free-energy formula */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if ( fabs(x) < 1e-8 ) {
      fe += -x*x*x/3 * ki2[i];
    } else {
      fe += (log(1 - x) + x + x*x/2) * ki2[i];
    }
    x = tk[i];
    //*ddmu -= x*x*x*ki2[i];
  }
  return fe * .5;
}



static void integ(int npt, double rmax, double rho)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk, *bphi, x;
  xdouble fe1, fe2, mu1, mu2, compr1, compr2, ddmu1, ddmu2;
  xdouble *fr, *dfr, *cr, *tr, *ck, *tk;
  xdouble *arr, *ri, *ki, *ri2, *ki2;
  int i;
  FFTWPFX(plan) plan = NULL;

  dr = rmax / npt;
  dk = PI / (dr * npt);

  xnew(arr, npt);
  xnew(ri, npt);
  xnew(ki, npt);
  xnew(ri2, npt);
  xnew(ki2, npt);
  plan = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);

  for ( i = 0; i < npt; i++ ) {
    ri[i] = dr * (i * 2 + 1) / 2;
    ki[i] = dk * (i * 2 + 1) / 2;
  }

  facr2k = PI*2 * dr;
  fack2r = pow(PI*2, -2) * dk;

  surfr = PI*4;
  surfk = PI*4 / pow(PI*2, 3);
  for ( i = 0; i < npt; i++ ) {
    ri2[i] = surfr * ri[i] * ri[i] * dr;
    ki2[i] = surfk * ki[i] * ki[i] * dk;
  }

  xnew(bphi, npt);
  xnew(fr, npt);
  xnew(dfr, npt);
  xnew(cr, npt);
  xnew(tr, npt);
  xnew(ck, npt);
  xnew(tk, npt);

  /* solvent-solvent interaction */
  for ( i = 0; i < npt; i++ ) {
    bphi[i] = beta * pot(ri[i], sigvv, 1, &dfr[i]);
    x = exp(-bphi[i]);
    fr[i] = x - 1;
    dfr[i] *= beta * x;
  }

  iter(npt, rho, ri, ki, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itermax);
  output(npt, ri, cr, tr, bphi, "vv1.dat");
  fe1 = getfe_hnc(npt, rho, cr, tr, ri2, ck, tk, ki2,
      &mu1, &compr1, &ddmu1);

  iter(npt, rho + delta, ri, ki, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itermax);
  output(npt, ri, cr, tr, bphi, "vv2.dat");
  fe2 = getfe_hnc(npt, rho + delta, cr, tr, ri2, ck, tk, ki2,
      &mu2, &compr2, &ddmu2);

  printf("rho %5.3f, T %6.3f: F1 %9.6f, F2 %9.6f, "
         "mu1 %9.6f, mu2 %9.6f, dF %9.6f; "
         "compr1 %9.6f, compr2 %9.6f, dmu(diff) %9.6f\n",
      (double) rho, (double) (1/beta),
      (double) fe1, (double) fe2, (double) mu1, (double) mu2,
      (double) ((fe2 - fe1)/delta),
      (double) compr1, (double) compr2,
      (double) ((mu2 - mu1)/delta) );
  printf("ddmu1 %9.6f, ddmu2 %9.6f, ddmu(diff) %9.6f\n",
      (double) ddmu1, (double) ddmu2,
      (double) ((compr2 - compr1)/delta) );

  free(bphi);
  free(fr);
  free(dfr);
  free(cr);
  free(tr);
  free(ck);
  free(tk);

  free(arr);
  free(ri);
  free(ki);
  free(ri2);
}



int main(int argc, char **argv)
{
  double den;

  doargs(argc, argv);

  for (den = 0.05; den <= rho; den += 0.05) {
    integ(numpt, rmax, den);
  }
#ifndef NOFFTW
  FFTWPFX(cleanup)();
#endif
  return 0;
}

