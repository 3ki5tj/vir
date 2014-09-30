/* free energy related quantities from the HNC closure */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



int dim = D;
int numpt = 32768;
int ffttype = 1;
xdouble rmax = (xdouble) 81.92L;
xdouble T = (xdouble) 3;
xdouble beta;
xdouble rhomax = (xdouble) 0.8L;
xdouble rhodel = (xdouble) 0.05L;
int itermax = 10000;
xdouble tol = (xdouble) 1e-8L;
xdouble damp = 1;
xdouble delta = (xdouble) 0.0001L;
int savecr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing free energy related quantities from the HNC closure";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rhomax, "maximal rho");
  argopt_add(ao, "--drho", "%" XDBLSCNF "f", &rhodel, "step size of rho");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--damp", "%" XDBLSCNF "f", &damp, "damping factor of solving integral equations");
  argopt_add(ao, "--del", "%" XDBLSCNF "f", &delta, "delta for numeric differentiation");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  printf("rmax %f, rho %g, T %f\n",
      (double) rmax, (double) rhomax, (double) T);
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
      x = (fr[i] + 1) * EXP(tr[i]) - tr[i] - 1;
      if ((err = FABS(cr[i] - x)) > errmax) errmax = err;
      cr[i] += damp * (x - cr[i]);
    }
    if ( errmax < tol ) break;
  }
  //printf("iter %d errmax %g\n", it, (double) errmax);
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
  xdouble fe1, fe2, mu1, mu2, compr1, compr2, ddmu1, ddmu2;
  xdouble P1, P2, P1a, P2a;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  sphr_t *sphr;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);

  xnew(bphi, npt);
  xnew(fr, npt);
  xnew(rdfr, npt);
  xnew(cr, npt);
  xnew(tr, npt);
  xnew(ck, npt);
  xnew(tk, npt);

  mkfr(npt, beta, bphi, fr, rdfr, sphr->ri, sphr->dm, 0, 0, 1);
  COPY1DARR(cr, fr, npt);

  for ( rho = rhodel; rho <= rhomax + 1e-8; rho += rhodel ) {
    iter(sphr, rho - delta, cr, tr, ck, tk, fr, itermax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, bphi, "vv1.dat");
    fe1 = getfe_hnc(npt, rho - delta, cr, tr, sphr->rDm1, ck, tk, sphr->kDm1,
        rdfr, &mu1, &compr1, &ddmu1, &P1, &P1a);

    iter(sphr, rho, cr, tr, ck, tk, fr, itermax);
    if ( savecr )
      output(npt, sphr->ri, cr, tr, bphi, "vv2.dat");
    fe2 = getfe_hnc(npt, rho, cr, tr, sphr->rDm1, ck, tk, sphr->kDm1,
        rdfr, &mu2, &compr2, &ddmu2, &P2, &P2a);

    printf("rho %5.3f, T %6.3f: F1 %9.6f, F2 %9.6f, "
           "mu1 %9.6f, mu2 %9.6f, dF %9.6f; "
           "compr1 %9.6f, compr2 %9.6f, dmu(diff) %9.6f\n",
        (double) rho, (double) (1/beta),
        (double) fe1, (double) fe2, (double) mu1, (double) mu2,
        (double) ((fe2 - fe1)/delta),
        (double) compr1, (double) compr2,
        (double) ((mu2 - mu1)/delta) );
    printf("ddmu1 %9.6f, ddmu2 %9.6f, ddmu(diff) %9.6f, P %9.6f %9.6f\n",
        (double) ddmu1, (double) ddmu2,
        (double) ((compr2 - compr1)/delta),
        (double) P2, (double) P2a);
  }

  sphr_close(sphr);
  free(bphi);
  free(fr);
  free(rdfr);
  free(cr);
  free(tr);
  free(ck);
  free(tk);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  integ(numpt, rmax);
  return 0;
}

