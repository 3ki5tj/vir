/* integral equation */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_MDIIS, SOLVER_COUNT };



int dim = D;
int numpt = 8192;
int ffttype = 1;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 1.35;
xdouble beta;
xdouble rho = (xdouble) 0.279L;
int ietype = IETYPE_PY;
int solver = SOLVER_LMV;
int itmax = 10000;
xdouble tol = (xdouble) 1e-8L;
xdouble damp = 1;
const xdouble errmax = 1000;
int Mpt = 20;
int nbases = 5;
char *fncrout = NULL;
char *fncrinp = NULL;
int verbose = 0;
xdouble params = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "solving an integral equation";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rho, "maximal rho");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-N", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-s", "%" XDBLSCNF "f", &params, "parameter of the integral equation");
  argopt_add(ao, "--damp", "%" XDBLSCNF "f", &damp, "damping factor of solving integral equations");
  argopt_add(ao, "--itmax", "%d", &itmax, "maximal number of iterations");
  argopt_add(ao, "--tol", "%" XDBLSCNF "f", &tol, "error tolerance");
  argopt_add(ao, "-S", "%d", &solver, "integral equation solver");
  argopt_add(ao, "-M", "%d", &Mpt, "number of points for the LMV solver");
  argopt_add(ao, "-B", "%d", &nbases, "number of bases for the MDIIS solver");
  argopt_add(ao, "--cr", NULL, &fncrout, "output file for c(r)");
  argopt_add(ao, "--inp", NULL, &fncrinp, "input file for c(r)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_add(ao, "--verbose", "%d", &verbose, "set verbose level");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  fprintf(stderr, "rmax %f, rho %g, T %f\n",
      (double) rmax, (double) rho, (double) T);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute the new indirect correlation function
 * *dcr = partial cr / partial tr  */
static xdouble getcr(xdouble tr, xdouble fr, int ietype, xdouble s,
    xdouble *dcr, xdouble *w, xdouble *dw)
{
  xdouble u, z;

  if ( ietype == IETYPE_PY ) { /* PY */
    if ( dcr ) *dcr = fr;
    return fr * (1 + tr);
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = EXP(tr);
    if ( dcr ) *dcr = (1 + fr) * z;
    if ( w ) *w = z - 1 - tr;
    if ( dw ) *dw = z - 1;
    return (1 + fr) * z - 1 - tr;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = EXP(tr);
    if ( dcr ) *dcr = (1 + fr) * (1 + s * (z - 1)) - 1;
    if (w) *w = z - 1 - tr;
    if (dw) *dw = z - 1;
    return fr * (1 + tr) + (1 + fr) * s * (z - 1 - tr);
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    u = 1 + s * tr;
    if (u < XDBL_MIN) u = XDBL_MIN;
    z = POW(u, 1./s);
    if ( dcr ) *dcr = (1 + fr) * z/u - 1;
    if ( w ) *w = (-LOG(u) + 1 - 1/u) * z / (s * s);
    if ( dw ) *dw = (-LOG(u) + (1 - s) * (1 - 1/u)) * z / u / (s * s);
    return (1 + fr) * z - 1 - tr;
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    if ( dcr ) *dcr = 1 + s * tr;
    if (w) *w = .5 * tr * tr;
    if (dw) *dw = tr;
    return fr * (1 + tr) + (1 + fr) * s * .5 * tr * tr;
  } else {
    fprintf(stderr, "unknown closure %d\n", ietype);
    exit(1);
  }
}



/* compute the new dcr = d c(r) / d(rho)
 * *d_dcr = partial dcr / partial dtr  */
static xdouble getdcr(xdouble dtr, xdouble tr, xdouble fr,
    int ietype, xdouble s, xdouble *d_dcr)
{
  xdouble z;

  if ( ietype == IETYPE_PY ) { /* PY */
    z = fr;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = (1 + fr) * EXP(tr) - 1;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = (1 + fr) * s * (EXP(tr) - 1) + fr;
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    z = (1 + fr) * POW(1 + s * tr, 1./s - 1) - 1;
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    z = (1 + fr) * s * tr + fr;
  } else {
    fprintf(stderr, "unknown closure %d\n", ietype);
    exit(1);
  }
  if ( d_dcr ) *d_dcr = z;
  return z * dtr;
}



/* do a step of direct (Picard) iteration
 * compute the residue vector `res' if needed */
static xdouble step_picard(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *fr, xdouble *res, xdouble dmp)
{
  int i, npt = sphr->npt;
  xdouble x, dx, err;

  sphr_r2k(sphr, cr, ck);
  for ( i = 0; i < npt; i++ )
    tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
  sphr_k2r(sphr, tk, tr);
  for ( err = 0, i = 0; i < npt; i++ ) {
    x = getcr(tr[i], fr[i], ietype, params, NULL, NULL, NULL);
    dx = x - cr[i];
    if ( res != NULL ) res[i] = dx;
    if ( FABS(dx) > err ) err = FABS(dx);
    cr[i] += dmp * dx;
  }
  return err;
}



/* do a step of direct (Picard) iteration for d/d(rho) functions
 * compute the residue vector `res' if needed */
static xdouble stepd_picard(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *fr, xdouble *res, xdouble dmp)
{
  int i, npt = sphr->npt;
  xdouble x, dx, hk, err;

  sphr_r2k(sphr, dcr, dck);
  for ( i = 0; i < npt; i++ ) {
    hk = ck[i] + tk[i];
    dtk[i] = hk * hk + rho * hk * (2 + rho * hk) * dck[i];
  }
  sphr_k2r(sphr, dtk, dtr);
  for ( err = 0, i = 0; i < npt; i++ ) {
    x = getdcr(dtr[i], tr[i], fr[i], ietype, params, NULL);
    dx = x - dcr[i];
    if ( res != NULL ) res[i] = dx;
    if ( FABS(dx) > err ) err = FABS(dx);
    dcr[i] += dmp * dx;
  }
  return err;
}



/* direct (Picard iteration) */
static xdouble iter_picard(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *fr, xdouble dmp, int itmax, xdouble tol)
{
  int it;
  xdouble err;

  for ( it = 0; it < itmax; it++ ) {
    err = step_picard(sphr, rho, cr, tr, ck, tk, fr, NULL, dmp);
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "iter_picard finished in %d steps, err %g\n", it, (double) err);
  return err;
}



/* solve d/d(rho) functions by direct iteration */
static void iterd_picard(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    xdouble *fr, xdouble dmp, int itmax, xdouble tol)
{
  int it;
  xdouble err;

  for ( it = 0; it < itmax; it++ ) {
    err = stepd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, fr, NULL, dmp);
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "iterd_picard finished in %d steps, err %g\n", it, (double) err);
}



/* advanced integral equation solvers */
const xdouble errinf = 1e300;
#include "ielmv.h"
#include "iemdiis.h"



/* load a previous solution for refinement */
static int loadcr(int npt, xdouble *ri, xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, const char *fn)
{
  int i;
  FILE *fp;
  char ln[40960];
  xdouble r, c, t, dc, dt;

  xfopen(fp, fn, "r", return -1);
  for ( i = 0; i < npt; i++ ) {
    /* get a data line */
    do {
      if ( fgets(ln, sizeof ln, fp) == NULL ) {
        fprintf(stderr, "%s corrupted for i %d\n", fn, i);
        goto EXIT;
      }
      strip(ln);
    } while ( ln[0] == '\0' || ln[0] == '#' );

    if ( 5 != sscanf(ln,
          "%" XDBLSCNF "f %" XDBLSCNF "f %" XDBLSCNF "f %" XDBLSCNF "f %" XDBLSCNF "f",
          &r, &c, &t, &dc, &dt) ) {
      fprintf(stderr, "%s no data on line %d\n%s\n", fn, i, ln);
      break;
    }
    if ( FABS(r - ri[i]) > 1e-6 ) {
      fprintf(stderr, "%s: %d, %g %g mismatch\n",
          fn, i, r, (double) ri[i]);
      break;
    }
    cr[i] = c;
    tr[i] = t;
    dcr[i] = dc;
    dtr[i] = dt;
  }
  fprintf(stderr, "%s loaded successfully\n", fn);
EXIT:
  fclose(fp);
  return i < npt;
}



static void savecr(int npt, xdouble *ri,
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



/* compute the beta pressure */
static xdouble getpres_py(int npt, xdouble rho,
    xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *rdfr,
    xdouble *ck,
    xdouble *ri2, xdouble *ki2,
    xdouble *ri, xdouble *ki,
    xdouble *compr, xdouble *dcompr,
    xdouble *presb, xdouble *comprb, xdouble *dcomprb,
    xdouble *presa)
{
  int i;
  xdouble pres, x;
  int once = 0;

  /* beta P = - rho^2/2 c(k = 0) + rho (1 + c(r = 0))/2
   *          + Int log(1 - rho c(k)) dk/(2pi)^3 } */
  pres = rho + ( -1 + get_zerosep(cr, ri) - rho * get_zerosep(ck, ki) ) * rho/2;
  *presa = rho -rho * rho * get_zerosep(ck, ki);
  *compr = 1;
  *dcompr = 0;
  *dcomprb = 0;
  for ( i = 0; i < npt; i++ ) {
    /* beta P = + rho^2/2 Int (h(r) - 1) c(r) dr
     *          + Int log(1 - rho c(k)) dk/(2pi)^3 } + rho c(r = 0)
     *        = + rho^2/2 Int (t(r) - 1) c(r) dr
     *          + Int {log(1 - rho c(k)) + rho c(k) + rho c(k)^2/2} dk/(2pi)^3 */
    /* compr = d(beta P) / drho = 1 - rho Int c(r) dr */
    *compr -= rho * cr[i] * ri2[i];
    /* exact value */
    *dcompr -= (cr[i] + rho * dcr[i]) * ri2[i];
    /* dcompr = d^2(beta P) / d(rho)^2 = -Int [c(r) + t(r) h(r)] dr
     * approximate result */
    *dcomprb -= (cr[i] + tr[i] * (cr[i] + tr[i])) * ri2[i];
    *presa -= rho * rho * (tr[i] + 1) * (tr[i] + 1) * rdfr[i] * ri2[i] / (2*dim);
  }
  /* the k-space part of the pressure formula */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if ( FABS(x) < 1e-8 ) {
      pres += (-x + x*x/2 - x*x*x/3) * ki2[i];
    } else if ( x < 1 ) {
      pres += LOG(1 - x) * ki2[i];
    } else {
      if ( !once )
        fprintf(stderr, "error for Pc: k %g, rho c %g\n", (double) ki[i], (double) x);
      once = 1;
    }
  }

  /* below are from the virial route */
  *presb = rho;
  *comprb = 1;
  for ( i = 0; i < npt; i++ ) {
    *presb += rdfr[i] * (rho*rho/6) * (1 + tr[i]) * ri2[i];
    *comprb += rdfr[i] * (rho/3 * (1 + tr[i]) + (rho*rho/6) * dtr[i]) * ri2[i];
  }
  return pres;
}



#if 0
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
  }
  fe *= 0.5;
  *pres = rho + rho * (*mu) - fe;
  *pres1 = rho + rho * rho * vir / (2 * dim);
  return fe;
}
#endif



static void integ(int npt, xdouble rmax)
{
  xdouble pres, presa, presb;
  xdouble compr, dcompr, comprb, dcomprb;
  xdouble *bphi, *der, *fr, *rdfr, *fk;
  xdouble *cr, *ck, *tr, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
  sphr_t *sphr;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);

  xnew(bphi, npt);
  xnew(der, npt);
  xnew(fr, npt);
  xnew(rdfr, npt);
  xnew(cr, npt);
  xnew(tr, npt);
  xnew(fk, npt);
  xnew(ck, npt);
  xnew(tk, npt);
  xnew(dcr, npt);
  xnew(dtr, npt);
  xnew(dck, npt);
  xnew(dtk, npt);

  mkfr(npt, beta, bphi, fr, rdfr, sphr->ri, sphr->dm, 0, 0, 1);
  sphr_r2k(sphr, fr, fk);

  if ( fncrinp != NULL ) {
    if ( loadcr(npt, sphr->ri, cr, tr, dcr, dtr, fncrinp) != 0 )
      goto EXIT;
    savecr(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, "initcr.dat");
  }

  if ( solver == SOLVER_LMV ) {
    iter_lmv(sphr, rho, cr, ck, tr, tk, fr, der, Mpt, damp, itmax, tol);
    iterd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, fr, damp, itmax, tol);
  } else if ( solver == SOLVER_MDIIS ) {
    iter_mdiis(sphr, rho, cr, tr, ck, tk, fr, nbases, damp, itmax, tol);
    iterd_mdiis(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, fr, nbases, damp, itmax, tol);
  } else {
    iter_picard(sphr, rho, cr, tr, ck, tk, fr, damp, itmax, tol);
    iterd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, fr, damp, itmax, tol);
  }
  if ( fncrout != NULL )
    savecr(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, fncrout);
  pres = getpres_py(npt, rho, cr, tr, dcr, dtr, rdfr,
      ck, sphr->rDm1, sphr->kDm1, sphr->ri, sphr->ki,
      &compr, &dcompr, &presb, &comprb, &dcomprb, &presa);

  printf("rho %5.3f, T %6.3f: Pc %9.5f/%9.5f, dPc %9.5f, ddPc %9.5f, "
      "Pv %9.5f, dPv %9.5f, ddPv %9.5f\n",
      (double) rho, (double) (1/beta),
      (double) pres, (double) presa, (double) compr, (double) dcompr,
      (double) presb, (double) comprb, (double) dcomprb);
EXIT:
  sphr_close(sphr);
  free(bphi);
  free(der);
  free(fr);
  free(rdfr);
  free(cr);
  free(tr);
  free(fk);
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

