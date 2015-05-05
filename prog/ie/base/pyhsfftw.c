/* pressure related quantities from the PY closure */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"
#include "ieutil.h"



int dim = D;
int numpt = 32768;
int ffttype = 1;
xdouble rmax = (xdouble) 81.92L;
xdouble rhomax = (xdouble) 0.8L;
xdouble rhodel = (xdouble) 0.05L;
int itmax = 10000;
xdouble tol = (xdouble) 1e-8L;
xdouble damp = 1;
xdouble delta = (xdouble) 0.0001L;
char *fnout = "PYhs.dat";
int savecr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing pressure and its derivatives from the PY closure";
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rhomax, "maxial rho");
  argopt_add(ao, "--drho", "%" XDBLSCNF "f", &rhodel, "step size of rho");
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--damp", "%" XDBLSCNF "f", &damp, "damping factor of solving integral equations");
  argopt_add(ao, "--del", "%" XDBLSCNF "f", &delta, "delta for numeric differentiation");
  argopt_add(ao, "-O", "%b", &savecr, "save cr and tr");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  printf("rmax %f, rho %g (%g)\n",
      (double) rmax, (double) rhomax, (double) delta);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



static void iter(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *fr, int itmax)
{
  int i, it, npt = sphr->npt;
  xdouble x, err, errmax;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, cr, ck);
    for ( i = 0; i < npt; i++ ) {
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    }
    sphr_k2r(sphr, tk, tr);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = fr[i] * (1 + tr[i]);
      if ((err = FABS(cr[i] - x)) > errmax) errmax = err;
      cr[i] += damp * (x - cr[i]);
    }
    if ( errmax < tol ) break;
  }
  //printf("iter %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(rho) functions */
static void iterd(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *ck, const xdouble *tk,
    xdouble *fr, int itmax)
{
  int i, it, npt = sphr->npt;
  xdouble x, hk, err, errmax;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, dcr, dck);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dtk[i] = hk*hk + rho*hk*(2 + rho*hk)*dck[i];
    }
    sphr_k2r(sphr, dtk, dtr);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      x = fr[i] * dtr[i];
      if ((err = FABS(dcr[i] - x)) > errmax) errmax = err;
      dcr[i] += damp * (x - dcr[i]);
    }
    //printf("round %d, errmax %g\n", it, errmax); getchar();
    if ( errmax < tol ) break;
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



/* compute the pressure */
static xdouble getpres_py(int npt, xdouble rho,
    xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *ri2,
    xdouble *ck, xdouble *ki2, xdouble *ri, xdouble *ki,
    xdouble B2, int dm,
    xdouble *compr, xdouble *dcompr,
    xdouble *presb, xdouble *comprb, xdouble *dcomprb,
    xdouble *pres1, xdouble *pres2)
{
  int i;
  xdouble pres, x;

  /* beta P = - rho^2/2 c(k = 0) - rho t(r = 0)/2
   *          + Int log(1 - rho c(k)) dk/(2pi)^3 } */
  pres = -get_zerosep(tr, ri)*rho/2 - get_zerosep(ck, ki)*rho*rho/2;
  /* beta P = rho - rho^2 Int c(r) dr - rho^2 Int [1 + t(r)]^2/(2D) r f'(r) \, dr */
  *pres1 = rho - rho * rho * B2 * (
      (tr[dm-1]+1)*(tr[dm-1]+1) + (tr[dm]+1)*(tr[dm]+1) ) * .5;
  /* beta P = rho - rho^2 c(k = 0) - rho^2 c^2(1) B2 */
  *pres2 = rho - rho * rho * get_zerosep(ck, ki)
           -rho * rho * B2 *
           ( 1.5*cr[dm-1]*cr[dm-1] - 0.5*cr[dm-2]*cr[dm-2] );
  *compr = 0;
  *dcompr = 0;
  *dcomprb = 0;
  for ( i = 0; i < npt; i++ ) {
    /* beta P = + rho^2/2 Int (h(r) - 1) c(r) dr
     *          + Int log(1 - rho c(k)) dk/(2pi)^3 } + rho c(r = 0)
     *        = + rho^2/2 Int (t(r) - 1) c(r) dr
     *          + Int {log(1 - rho c(k)) + rho c(k) + rho c(k)^2/2} dk/(2pi)^3 */
    /* compr = d(beta P) / drho = -rho Int c(r) dr */
    *compr -= cr[i] * ri2[i];
    /* exact value */
    *dcompr -= (cr[i] + rho * dcr[i]) * ri2[i];
    /* dcompr = d^2(beta P) / d(rho)^2 = -Int [c(r) + t(r) h(r)] dr
     * approximate result */
    *dcomprb -= (cr[i] + tr[i] * (cr[i] + tr[i])) * ri2[i];
    *pres1 -= rho * rho * cr[i] * ri2[i];
  }
  *compr *= rho;
  /* the k-space part of the pressure formula */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if ( FABS(x) < 1e-8 ) {
      pres += x*(-1 + x*(1./2 - x/3)) * ki2[i];
    } else {
      pres += LOG(1 - x) * ki2[i];
    }
  }

  /* below are from the pressure route */
  *presb = rho + B2 * (rho*rho) * (1 + tr[dm]);
  *comprb = B2 * rho * (2 + 2*tr[dm] + rho * dtr[dm]);
  return pres;
}



static void integ(int npt, xdouble rmax)
{
  xdouble rho, *bphi;
  xdouble pres1, pres2, compr1, compr2, dcompr1, dcompr2;
  xdouble pres1b, pres2b, compr1b, compr2b, dcompr1b, dcompr2b;
  xdouble pres1u, pres1v, pres2u, pres2v;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
  sphr_t *sphr;
  FILE *fp;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);

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

  mkfr(npt, 1, bphi, NULL, NULL, fr, rdfr, sphr->ri, sphr->dm, 0, 0, 0);
  COPY1DARR(cr, fr, npt);
  COPY1DARR(dcr, fr, npt);

  xfopen(fp, fnout, "w", exit(1));

  for ( rho = rhodel; rho <= rhomax + 1e-8; rho += rhodel ) {
    /* 1. at the density of rho - delta */
    /* 1A. compute c(r), t(r) ... */
    iter(sphr, rho - delta, cr, tr, ck, tk, fr, itmax);
    /* 1B. compute dc(r), dt(r) ... */
    iterd(sphr, rho - delta, dcr, dtr, dck, dtk, ck, tk,
        fr, itmax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, "vv1.dat");
    pres1 = getpres_py(npt, rho - delta, cr, tr, dcr, dtr,
        sphr->rDm1, ck, sphr->kDm1, sphr->ri, sphr->ki,
        sphr->B2hs, sphr->dm, &compr1, &dcompr1, &pres1b, &compr1b, &dcompr1b,
        &pres1u, &pres1v);

    /* at the density of rho */
    /* 2. at the density of rho */
    /* 2A. compute c(r), t(r) ... */
    iter(sphr, rho, cr, tr, ck, tk, fr, itmax);
    /* 2B. compute dc(r), dt(r) ... */
    iterd(sphr, rho, dcr, dtr, dck, dtk, ck, tk, fr, itmax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, "vv2.dat");
    pres2 = getpres_py(npt, rho, cr, tr, dcr, dtr,
        sphr->rDm1, ck, sphr->kDm1, sphr->ri, sphr->ki,
        sphr->B2hs, sphr->dm, &compr2, &dcompr2, &pres2b, &compr2b, &dcompr2b,
        &pres2u, &pres2v);

    printf("rho %5.3f, Pc- %8.5f, Pc %8.5f, "
           "dPc- %8.5f, dPc %8.5f, dP %8.5f; "
           "ddPc- %8.4f, ddPc %8.4f, ddP(diff) %8.4f | Pc %8.5f, %8.5f\n",
        (double) rho,
        (double) pres1, (double) pres2,
        (double) compr1, (double) compr2,
        (double) ((pres2 - pres1)/delta),
        (double) dcompr1, (double) dcompr2,
        (double) ((compr2 - compr1)/delta),
        (double) pres2u, (double) pres2v);

    printf("rho %5.3f, Pv- %8.5f, Pv %8.5f, "
           "dPv- %8.5f, dPv %8.5f, dP %8.5f; "
           "ddP-  %8.4f, ddP  %8.4f, ddP(diff) %8.4f\n",
        (double) rho,
        (double) pres1b, (double) pres2b,
        (double) compr1b, (double) compr2b,
        (double) ((pres2b - pres1b)/delta),
        (double) dcompr1b, (double) dcompr2b,
        (double) ((compr2b - compr1b)/delta) );

    printf("\n");
    fprintf(fp, "%g %g %g %g %g\n", (double) rho, (double) pres2,
        (double) pres2u, (double) pres2v, (double) pres2b);
  }

  fclose(fp);
  sphr_close(sphr);
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
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  integ(numpt, rmax);
  return 0;
}

