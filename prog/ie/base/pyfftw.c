/* pressure related quantities from the PY closure */
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
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



int numpt = 8192;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 10;
xdouble beta;
xdouble rho = (xdouble) 0.8L;
int itmax = 10000;
xdouble tol = (xdouble) 1e-12L;
xdouble delta = (xdouble) 0.001L;
int savecr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing pressure and its derivatives from the HNC closure";
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "-r", "%" XDBLSCNF "f", &rho, "rho");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &delta, "delta lambda");
  argopt_add(ao, "-O", "%b", &savecr, "save cr and tr");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  printf("rmax %f, rho %g (%g), T %f\n",
      (double) rmax, (double) rho, (double) delta, (double) T);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* get the value at zero separation */
__inline static xdouble get_zerosep(xdouble *y, xdouble *x)
{
  return (y[0]*x[1] - y[1]*x[0])/(x[1] - x[0]);
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
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, err, errmax;

  for ( i = 0; i < npt; i++ ) cr[i] = fr[i];
  for ( it = 0; it < itmax; it++ ) {
    sphr(npt, cr, ck, facr2k, plan, arr, ri, ki);
    for ( i = 0; i < npt; i++ ) {
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    }
    sphr(npt, tk, tr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = fr[i] * (1 + tr[i]);
      if ((err = FABS(cr[i] - x)) > errmax) errmax = err;
      cr[i] = x;
    }
    if ( errmax < 1e-8 ) break;
  }
  //printf("iter %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(rho) functions */
static void iterd(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *cr, const xdouble *tr, const xdouble *ck, const xdouble *tk,
    xdouble *fr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, hk, err, errmax;

  for ( i = 0; i < npt; i++ ) dcr[i] = fr[i];
  for ( it = 0; it < itmax; it++ ) {
    sphr(npt, dcr, dck, facr2k, plan, arr, ri, ki);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dtk[i] = hk*hk + rho*hk*(2 + rho*hk)*dck[i];
    }
    sphr(npt, dtk, dtr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = fr[i] * dtr[i];
      if ((err = FABS(dcr[i] - x)) > errmax) errmax = err;
      dcr[i] = x;
    }
    //printf("round %d, errmax %g\n", it, errmax); getchar();
    if ( errmax < 1e-7 ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



static void output(int npt, xdouble *ri,
    xdouble *cr, xdouble *tr, xdouble *dcr, xdouble *dtr,
    xdouble *fr, xdouble *bphi, char *fn)
{
  int i;
  FILE *fp;

  xfopen(fp, fn, "w", return);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%8.6f %14.8f %14.8f %14.8f %14.8f %14.8f %g\n",
        (double) ri[i],
        (double) cr[i], (double) tr[i],
        (double) dcr[i], (double) dtr[i],
        (double) fr[i], (double) bphi[i]);
  }
  fclose(fp);
}



/* compute the excess pressure */
static xdouble getpres_py(int npt, xdouble rho,
    xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *rdfr, xdouble *ri2,
    xdouble *ck, xdouble *ki2, xdouble *ri, xdouble *ki,
    xdouble *compr, xdouble *dcompr,
    xdouble *presb, xdouble *comprb, xdouble *dcomprb)
{
  int i;
  xdouble pres, x;

  //pres = rho;
  /* beta P = - rho^2/2 c(k = 0) + rho (1 + c(r = 0))/2
   *          + Int log(1 - rho c(k)) dk/(2pi)^3 } */
  pres = (-1 + get_zerosep(cr, ri) - rho * get_zerosep(ck, ki) ) * rho/2;
  *compr = 0;
  *dcompr = 0;
  *dcomprb = 0;
  for ( i = 0; i < npt; i++ ) {
    /* beta P = + rho^2/2 Int (h(r) - 1) c(r) dr
     *          + Int log(1 - rho c(k)) dk/(2pi)^3 } + rho c(r = 0)
     *        = + rho^2/2 Int (t(r) - 1) c(r) dr
     *          + Int {log(1 - rho c(k)) + rho c(k) + rho c(k)^2/2} dk/(2pi)^3 */
    //pres += 0.5 * cr[i] * (tr[i] - 1) * ri2[i];
    /* compr = d(beta P) / drho = -rho Int c(r) dr */
    *compr -= cr[i] * ri2[i];
    /* exact value */
    *dcompr -= (cr[i] + rho * dcr[i]) * ri2[i];
    /* dcompr = d^2(beta P) / d(rho)^2 = -Int [c(r) + t(r) h(r)] dr
     * approximate result */
    *dcomprb -= (cr[i] + tr[i] * (cr[i] + tr[i])) * ri2[i];
  }
  *compr *= rho;
  /* the k-space part of the pressure formula */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    //if ( FABS(x) < 1e-8 ) {
    //  pres += - x*x*x/3 * ki2[i];
    //} else {
    //  pres += (LOG(1 - x) + x - x*x/2)* ki2[i];
    //}
    if ( FABS(x) < 1e-8 ) {
      pres += (-x + x*x/2 - x*x*x/3) * ki2[i];
    } else {
      pres += LOG(1 - x) * ki2[i];
    }
  }

  /* below are from the pressure route */
  *presb = 0;
  *comprb = 0;
  for ( i = 0; i < npt; i++ ) {
    *presb += rdfr[i] * (rho*rho/6) * (1 + tr[i]) * ri2[i];
    *comprb += rdfr[i] * (rho/3 * (1 + tr[i]) + (rho*rho/6) * dtr[i]) * ri2[i];
  }
  return pres;
}



static void integ(int npt, xdouble rmax, xdouble rho)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk, *bphi, x;
  xdouble pres1, pres2, compr1, compr2, dcompr1, dcompr2;
  xdouble pres1b, pres2b, compr1b, compr2b, dcompr1b, dcompr2b;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
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
  fack2r = pow_si(PI*2, -2) * dk;

  surfr = PI*4;
  surfk = surfr * pow_si(PI*2, -3);
  for ( i = 0; i < npt; i++ ) {
    ri2[i] = surfr * ri[i] * ri[i] * dr;
    ki2[i] = surfk * ki[i] * ki[i] * dk;
  }

  xnew(bphi, npt);
  xnew(fr, npt);
  xnew(rdfr, npt);
  xnew(cr, npt);
  xnew(tr, npt);
  xnew(ck, npt);
  xnew(tk, npt);
  xnew(dcr, npt);
  xnew(dtr, npt);
  xnew(dck, npt);
  xnew(dtk, npt);

  /* solvent-solvent interaction */
  for ( i = 0; i < npt; i++ ) {
    bphi[i] = beta * pot(ri[i], 1, 1, &rdfr[i]);
    x = EXP(-bphi[i]);
    fr[i] = x - 1;
    rdfr[i] *= beta * x;
  }

  iter(npt, rho - delta, ri, ki, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itmax);
  iterd(npt, rho, ri, ki, dcr, dtr, dck, dtk, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itmax);
  if ( savecr )
    output(npt, ri, cr, tr, dcr, dtr, fr, bphi, "vv1.dat");
  pres1 = getpres_py(npt, rho - delta, cr, tr, dcr, dtr, rdfr, ri2, ck, ki2, ri, ki,
      &compr1, &dcompr1, &pres1b, &compr1b, &dcompr1b);

  iter(npt, rho, ri, ki, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itmax);
  iterd(npt, rho, ri, ki, dcr, dtr, dck, dtk, cr, tr, ck, tk,
      fr, facr2k, fack2r, plan, arr, itmax);
  if ( savecr )
    output(npt, ri, cr, tr, dcr, dtr, fr, bphi, "vv2.dat");
  pres2 = getpres_py(npt, rho, cr, tr, dcr, dtr, rdfr, ri2, ck, ki2, ri, ki,
      &compr2, &dcompr2, &pres2b, &compr2b, &dcompr2b);


  printf("rho %5.3f, T %6.3f: P1 %9.6f, P2 %9.6f, "
         "dP1 %9.6f, dP2 %9.6f, dP %9.6f; "
         "ddP1 %9.6f, ddP2 %9.6f, ddP(diff) %9.6f\n",
      (double) rho, (double) (1/beta),
      (double) pres1, (double) pres2, (double) compr1, (double) compr2,
      (double) ((pres2 - pres1)/delta),
      (double) dcompr1, (double) dcompr2,
      (double) ((compr2 - compr1)/delta) );

  printf("rho %5.3f, T %6.3f: P1 %9.6f, P2 %9.6f, "
         "dP1 %9.6f, dP2 %9.6f, dP %9.6f; "
         "ddP1 %9.6f, ddP2 %9.6f, ddP(diff) %9.6f\n\n",
      (double) rho, (double) (1/beta),
      (double) pres1b, (double) pres2b, (double) compr1b, (double) compr2b,
      (double) ((pres2b - pres1b)/delta),
      (double) dcompr1b, (double) dcompr2b,
      (double) ((compr2b - compr1b)/delta) );

  free(bphi);
  free(fr);
  free(rdfr);
  free(cr);
  free(tr);
  free(ck);
  free(tk);
  free(dcr);
  free(dtr);
  free(dck);
  free(dtk);

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

