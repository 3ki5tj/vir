/* integral equation */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



int dim = D;
int numpt = 8192;
int ffttype = 1;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 3;
xdouble beta;
xdouble rhomax = (xdouble) 0.8L;
xdouble rhodel = (xdouble) 0.05L;
int itmax = 10000;
xdouble tol = (xdouble) 1e-10L;
xdouble damp = 1;
xdouble delta = (xdouble) 0.001L;
int savecr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "solving an integral equation";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rhomax, "maximal rho");
  argopt_add(ao, "--drho", "%" XDBLSCNF "f", &rhodel, "step size of rho");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--damp", "%" XDBLSCNF "f", &damp, "damping factor of solving integral equations");
  argopt_add(ao, "--del", "%" XDBLSCNF "f", &delta, "delta for numeric differentiation");
  argopt_add(ao, "-O", "%b", &savecr, "save cr and tr");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  fprintf(stderr, "rmax %f, rho %g (%g), T %f\n",
      (double) rmax, (double) rhomax, (double) delta, (double) T);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



static double getcr(double tr, double fr, int ietype)
{
  if ( ietype == IETYPE_PY ) {
    return fr * (1 + tr);
  } else if ( ietype == IETYPE_HNC ) {
    return (1 + fr) * EXP( tr ) - 1 - tr;
  } else {
    fprintf(stderr, "unknown closure\n");
    exit(1);
  }
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
      x = getcr(tr[i], fr[i], ietype);
      if ((err = FABS(x - cr[i])) > errmax) errmax = err;
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
    xdouble *presb, xdouble *comprb, xdouble *dcomprb,
    xdouble *presa)
{
  int i;
  xdouble pres, x;

  //pres = rho;
  /* beta P = - rho^2/2 c(k = 0) + rho (1 + c(r = 0))/2
   *          + Int log(1 - rho c(k)) dk/(2pi)^3 } */
  pres = (-1 + get_zerosep(cr, ri) - rho * get_zerosep(ck, ki) ) * rho/2;
  *presa = -rho * rho * get_zerosep(ck, ki);
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
    *presa -= rho * rho * (tr[i] + 1) * (tr[i] + 1) * rdfr[i] * ri2[i] / (2*dim);
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

  /* below are from the virial route */
  *presb = 0;
  *comprb = 0;
  for ( i = 0; i < npt; i++ ) {
    *presb += rdfr[i] * (rho*rho/6) * (1 + tr[i]) * ri2[i];
    *comprb += rdfr[i] * (rho/3 * (1 + tr[i]) + (rho*rho/6) * dtr[i]) * ri2[i];
  }
  return pres;
}



static xdouble getfe_hnc(int npt, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ri2,
    xdouble *ck, xdouble *tk, xdouble *ki2,
    xdouble *rdfr,
    xdouble *mu, xdouble *compr, xdouble *ddmu,
    xdouble *pres, xdouble *pres1)
{
  int i;
  xdouble fe, x, vir = 0;

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
    *ddmu -= tr[i] * (tr[i] + cr[i]) * ri2[i];
    vir += EXP(tr[i]) * rdfr[i] * ri2[i];
  }
  fe *= rho * rho;
  *mu *= rho;
  *ddmu /= rho;
  /* the k-space part of the free-energy formula */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if ( FABS(x) < 1e-8 ) {
      fe += -x*x*x/3 * ki2[i];
    } else {
      fe += (LOG(1 - x) + x + x*x/2) * ki2[i];
    }
    x = tk[i];
    //*ddmu -= x*x*x*ki2[i];
  }
  fe *= 0.5;
  *pres = rho + rho * (*mu) - fe;
  *pres1 = rho + rho * rho * vir / (2 * dim);
  return fe;
}



static void integ(int npt, xdouble rmax)
{
  xdouble rho, *bphi;
  xdouble pres1, pres2, compr1, compr2, dcompr1, dcompr2;
  xdouble pres1b, pres2b, compr1b, compr2b, dcompr1b, dcompr2b;
  xdouble pres1a, pres2a;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
  sphr_t *sphr;

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

  mkfr(npt, beta, bphi, fr, rdfr, sphr->ri, sphr->dm, 0, 0, 1);

  for ( rho = rhodel; rho <= rhomax + 1e-8; rho += rhodel ) {
    iter(sphr, rho - delta, cr, tr, ck, tk, fr, itmax);
    iterd(sphr, rho - delta, dcr, dtr, dck, dtk, ck, tk, fr, itmax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, "vv1.dat");
    pres1 = getpres_py(npt, rho - delta, cr, tr, dcr, dtr, rdfr,
        sphr->rDm1, ck, sphr->kDm1, sphr->ri, sphr->ki,
        &compr1, &dcompr1, &pres1b, &compr1b, &dcompr1b,
        &pres1a);

    iter(sphr, rho, cr, tr, ck, tk, fr, itmax);
    iterd(sphr, rho, dcr, dtr, dck, dtk, ck, tk, fr, itmax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, "vv2.dat");
    pres2 = getpres_py(npt, rho, cr, tr, dcr, dtr, rdfr,
        sphr->rDm1, ck, sphr->kDm1, sphr->ri, sphr->ki,
        &compr2, &dcompr2, &pres2b, &compr2b, &dcompr2b,
        &pres2a);

    printf("rho %5.3f, T %6.3f: P1 %9.6f, P2 %9.6f, "
           "dP1 %9.6f, dP2 %9.6f, dP %9.6f; "
           "ddP1 %9.6f, ddP2 %9.6f, ddP(diff) %9.6f | Pc %8.5f\n",
        (double) rho, (double) (1/beta),
        (double) pres1, (double) pres2, (double) compr1, (double) compr2,
        (double) ((pres2 - pres1)/delta),
        (double) dcompr1, (double) dcompr2,
        (double) ((compr2 - compr1)/delta),
        (double) pres2a);

    printf("rho %5.3f, T %6.3f: Pv-%9.6f, Pv %9.6f, "
           "dPv-%9.6f, dPv %9.6f, dP %9.6f; "
           "ddPv-%9.6f, ddPv %9.6f, ddP(diff) %9.6f\n",
        (double) rho, (double) (1/beta),
        (double) pres1b, (double) pres2b, (double) compr1b, (double) compr2b,
        (double) ((pres2b - pres1b)/delta),
        (double) dcompr1b, (double) dcompr2b,
        (double) ((compr2b - compr1b)/delta));
  }

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

