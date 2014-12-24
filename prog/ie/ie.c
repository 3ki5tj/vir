/* integral equation
 * properties around a nonzero density */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_MDIIS, SOLVER_COUNT };
const char *solvers[] = {"Picard", "LMV", "MDIIS"};



int dim = D;
int numpt = 8192;
int ffttype = 1;
xdouble rmax = (xdouble) 40.96L;
xdouble T = (xdouble) 1.35;
xdouble beta;
xdouble rho = (xdouble) 0.279L;
int ljtype = LJ_FULL;
int ietype = IETYPE_PY;
int solver = SOLVER_MDIIS;
int itmax = 10000;
xdouble tol = (xdouble) 1e-8L;
xdouble damp = 1;
const xdouble errmax = 10000;
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
  argopt_addx(ao, "--LJ", "%list", &ljtype, "LJ type", ljtypes, LJ_COUNT);
  argopt_add(ao, "--LJlrs", "%" XDBLSCNF "f", &ljlrs, "scaling factor for the long-range part when using --LJ=split");
  argopt_addx(ao, "-S", "%list", &solver, "integral equation solver", solvers, SOLVER_COUNT);
  argopt_addx(ao, "-C", "%list", &ietype, "closure", ietype_names, IETYPE_COUNT);
  argopt_add(ao, "-M", "%d", &Mpt, "number of points for the LMV solver");
  argopt_add(ao, "-B", "%d", &nbases, "number of bases for the MDIIS solver");
  argopt_add(ao, "--cr", NULL, &fncrout, "output file for c(r)");
  argopt_add(ao, "--in", NULL, &fncrinp, "input file for c(r)");
  argopt_add(ao, "-v", "%+", &verbose, "be verbose");
  argopt_add(ao, "--verbose", "%d", &verbose, "set verbose level");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1/T;
  if ( Mpt <= 0 ) Mpt = (int) (2 * rmax);
  fprintf(stderr, "rmax %f, rho %g, T %f, closure %s, solver %s, "
      "M %d, nbases %d\n",
      (double) rmax, (double) rho, (double) T,
      solvers[solver], ietype_names[ietype], Mpt, nbases);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute tk from ck */
static void oz(xdouble *ck, xdouble *tk,
    xdouble *der, int npt, xdouble rho)
{
  int i;

  for ( i = 0; i < npt; i++ ) {
    xdouble y = rho * ck[i] / (1 - rho * ck[i]);
    tk[i] = ck[i] * y;
    if ( der != NULL ) der[i] = y * (2 + y);
  }
}



/* compute dtk from dck */
static void ozd(xdouble *dck, xdouble *dtk,
    const xdouble *ck, const xdouble *tk,
    xdouble *der, int npt, xdouble rho)
{
  int i;

  for ( i = 0; i < npt; i++ ) {
    xdouble hk = ck[i] + tk[i];
    xdouble K = 1 + rho * hk;
    xdouble dtdc = K*K - 1;
    dtk[i] = hk * hk + dtdc * dck[i];
    if ( der != NULL ) der[i] = dtdc;
  }
}



/* compute ddtk from ddck */
static void ozdd(xdouble *ddck, xdouble *ddtk,
    const xdouble *dck, const xdouble *dtk,
    const xdouble *ck, const xdouble *tk,
    xdouble *der, int npt, xdouble rho)
{
  int i;

  for ( i = 0; i < npt; i++ ) {
    xdouble hk = ck[i] + tk[i];
    xdouble K = 1 + rho * hk, dtdc = K*K - 1;
    xdouble dhk = dck[i] + dtk[i];
    ddtk[i] = K * dhk * ck[i] + K * dck[i] * (hk + rho * dhk) + dtdc * ddck[i];
    if ( der != NULL ) der[i] = dtdc;
  }
}



/* compute the cavity distribution function
 * *d_yr = partial yr / partial tr  */
static xdouble getyr(xdouble tr, xdouble bphilr, int ietype, xdouble s,
    xdouble *d_yr, xdouble *dd_yr)
{
  xdouble u, z, tr1;

  tr1 = tr - bphilr;
  if ( ietype == IETYPE_PY ) { /* PY */
    if ( d_yr ) *d_yr = 1;
    if ( dd_yr ) *dd_yr = 0;
    return  1 + tr1;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = EXP(tr1);
    if ( d_yr ) *d_yr = z;
    if ( dd_yr ) *dd_yr = z;
    return z;
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    if ( d_yr ) *d_yr = 1 + s * tr1;
    if ( dd_yr ) *dd_yr = s;
    return 1 + tr1 + s * .5 * tr1 * tr1;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = EXP(tr1);
    if ( d_yr ) *d_yr = 1 + s * (z - 1);
    if ( dd_yr ) *dd_yr = s * z;
    return 1 + tr1 + s * (z - 1 - tr1);
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    u = 1 + s * tr1;
    if (u < XDBL_MIN) u = XDBL_MIN;
    z = POW(u, 1./s);
    if ( d_yr ) *d_yr = z/u;
    if ( dd_yr ) *dd_yr = (1 - s) * z / (u * u);
    return z;
  } else {
    fprintf(stderr, "cannot do closure %s\n", ietype_names[ietype]);
    exit(1);
  }
}



/* compute the new indirect correlation function
 * *d_cr = partial cr / partial tr  */
static xdouble getcr(xdouble tr, xdouble bphilr,
    xdouble fr, int ietype, xdouble s,
    xdouble *d_cr, xdouble *w, xdouble *dw)
{
  xdouble u, z, tr1;

  tr1 = tr - bphilr;
  if ( ietype == IETYPE_PY ) { /* PY */
    if ( d_cr ) *d_cr = fr;
    return (1 + fr) * (1 + tr1) - 1 - tr;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = EXP(tr1);
    if ( d_cr ) *d_cr = (1 + fr) * z;
    if ( w ) *w = z - 1 - tr;
    if ( dw ) *dw = z - 1;
    return (1 + fr) * z - 1 - tr;
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    if ( d_cr ) *d_cr = (1 + fr) * (1 + s * tr1) - 1;
    if (w) *w = .5 * tr1 * tr1;
    if (dw) *dw = tr1;
    return (1 + fr) * (1 + tr1 + s * .5 * tr1 * tr1) - 1 - tr;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = EXP(tr1);
    if ( d_cr ) *d_cr = (1 + fr) * (1 + s * (z - 1)) - 1;
    if (w) *w = z - 1 - tr1;
    if (dw) *dw = z - 1;
    return (1 + fr) * (1 + tr1 + s * (z - 1 - tr1)) - 1 - tr;
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    u = 1 + s * tr1;
    if (u < XDBL_MIN) u = XDBL_MIN;
    z = POW(u, 1./s);
    if ( d_cr ) *d_cr = (1 + fr) * z/u - 1;
    if ( w ) *w = (-LOG(u) + 1 - 1/u) * z / (s * s);
    if ( dw ) *dw = (-LOG(u) + (1 - s) * (1 - 1/u)) * z / u / (s * s);
    return (1 + fr) * z - 1 - tr;
  } else {
    fprintf(stderr, "cannot do closure %s\n", ietype_names[ietype]);
    exit(1);
  }
}



/* compute the new dcr = d c(r) / d(rho)
 * *d_dcr = partial dcr / partial dtr  */
static xdouble getdcr(xdouble dtr, xdouble tr,
    xdouble bphilr, xdouble fr,
    int ietype, xdouble s, xdouble *d_dcr)
{
  xdouble z, tr1;

  tr1 = tr - bphilr;
  if ( ietype == IETYPE_PY ) { /* PY */
    z = fr;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = (1 + fr) * EXP(tr1) - 1;
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    z = (1 + fr) * s * tr1 + fr;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = (1 + fr) * s * (EXP(tr1) - 1) + fr;
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    z = (1 + fr) * POW(1 + s * tr1, 1./s - 1) - 1;
  } else {
    fprintf(stderr, "cannot do closure %s\n", ietype_names[ietype]);
    exit(1);
  }
  if ( d_dcr ) *d_dcr = z;
  return z * dtr;
}



/* compute the new ddcr = d^2 c(r) / d(rho)^2
 * *d_ddcr = partial dcr / partial ddtr  */
static xdouble getddcr(xdouble ddtr, xdouble dtr, xdouble tr,
    xdouble bphilr, xdouble fr,
    int ietype, xdouble s, xdouble *d_ddcr)
{
  xdouble xp, z = 0, w = 0, tr1;

  tr1 = tr - bphilr;
  if ( ietype == IETYPE_PY ) { /* PY */
    z = fr;
    w = 0;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    xp = EXP(tr1);
    z = (1 + fr) * xp - 1;
    w = (1 + fr) * xp;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    xp = EXP(tr1);
    z = (1 + fr) * s * (xp - 1) + fr;
    w = (1 + fr) * s * xp;
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    z = (1 + fr) * POW(1 + s * tr1, 1./s - 1) - 1;
    w = (1 + fr) * (1 - s) * POW(1 + s * tr1, 1./s - 2);
  } else if ( ietype == IETYPE_SQR ) { /* quadratic */
    z = (1 + fr) * s * tr1 + fr;
    w = (1 + fr) * s;
  } else {
    fprintf(stderr, "cannot do closure %s\n", ietype_names[ietype]);
    exit(1);
  }
  if ( d_ddcr ) *d_ddcr = z;
  return w * dtr * dtr + z * ddtr;
}



/* closure for correlation functions */
static xdouble closure(xdouble *cr, xdouble *tr,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble *der, int npt,
    int ietype, xdouble params, xdouble dmp)
{
  int i;
  xdouble err, x, dx, dc;

  for ( err = 0, i = 0; i < npt; i++ ) {
    x = getcr(tr[i], bphilr[i], fr[i], ietype, params, &dc, NULL, NULL);
    dx = x - cr[i];
    if ( der != NULL ) der[i] = dc;
    if ( res != NULL ) res[i] = dx;
    if ( FABS(dx) > err ) err = FABS(dx);
    cr[i] += dmp * dx;
  }
  return err;
}



/* closure for d/d(rho) functions */
static xdouble closured(xdouble *dcr, xdouble *dtr,
    const xdouble *tr,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble *der, int npt,
    int ietype, xdouble params, xdouble dmp)
{
  int i;
  xdouble err, x, dx, dc;

  for ( err = 0, i = 0; i < npt; i++ ) {
    x = getdcr(dtr[i], tr[i], bphilr[i], fr[i], ietype, params, &dc);
    dx = x - dcr[i];
    if ( der != NULL ) der[i] = dc;
    if ( res != NULL ) res[i] = dx;
    if ( FABS(dx) > err ) err = FABS(dx);
    dcr[i] += dmp * dx;
  }
  return err;
}



/* closure for d^2/d(rho)^2 functions */
static xdouble closuredd(xdouble *ddcr, xdouble *ddtr,
    const xdouble *dtr, const xdouble *tr,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble *der, int npt,
    int ietype, xdouble params, xdouble dmp)
{
  int i;
  xdouble err, x, dx, dc;

  for ( err = 0, i = 0; i < npt; i++ ) {
    x = getddcr(ddtr[i], dtr[i], tr[i], bphilr[i], fr[i], ietype, params, &dc);
    dx = x - ddcr[i];
    if ( der != NULL ) der[i] = dc;
    if ( res != NULL ) res[i] = dx;
    if ( FABS(dx) > err ) err = FABS(dx);
    ddcr[i] += dmp * dx;
  }
  return err;
}



/* do a step of direct (Picard) iteration
 * compute the residue vector `res' if needed */
static xdouble step_picard(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble dmp)
{
  sphr_r2k(sphr, cr, ck);
  oz(ck, tk, NULL, sphr->npt, rho);
  sphr_k2r(sphr, tk, tr);
  return closure(cr, tr, bphilr, fr, res, NULL, sphr->npt, ietype, params, dmp);
}



/* do a step of direct (Picard) iteration for d/d(rho) functions
 * compute the residue vector `res' if needed */
static xdouble stepd_picard(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble dmp)
{
  sphr_r2k(sphr, dcr, dck);
  ozd(dck, dtk, ck, tk, NULL, sphr->npt, rho);
  sphr_k2r(sphr, dtk, dtr);
  return closured(dcr, dtr, tr, bphilr, fr, res, NULL,
      sphr->npt, ietype, params, dmp);
}



/* do a step of direct (Picard) iteration for d^2/d(rho)^2 functions
 * compute the residue vector `res' if needed */
static xdouble stepdd_picard(sphr_t *sphr, xdouble rho,
    xdouble *ddcr, xdouble *ddtr, xdouble *ddck, xdouble *ddtk,
    const xdouble *dtr, const xdouble *dck, const xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble *res, xdouble dmp)
{
  sphr_r2k(sphr, ddcr, ddck);
  ozdd(ddck, ddtk, dck, dtk, ck, tk, NULL, sphr->npt, rho);
  sphr_k2r(sphr, ddtk, ddtr);
  return closuredd(ddcr, ddtr, dtr, tr, bphilr, fr, res, NULL,
      sphr->npt, ietype, params, dmp);
}



/* direct (Picard iteration) */
static xdouble iter_picard(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble dmp, int itmax, xdouble tol)
{
  int it;
  xdouble err = 0;

  for ( it = 0; it < itmax; it++ ) {
    err = step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, NULL, dmp);
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "iter_picard finished in %d steps, err %g\n", it, (double) err);
  return err;
}



/* solve d/d(rho) functions by direct iteration */
static xdouble iterd_picard(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble dmp, int itmax, xdouble tol)
{
  int it;
  xdouble err = 0;

  for ( it = 0; it < itmax; it++ ) {
    err = stepd_picard(sphr, rho, dcr, dtr, dck, dtk,
        tr, ck, tk, bphilr, fr, NULL, dmp);
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "iterd_picard finished in %d steps, err %g\n", it, (double) err);
  return err;
}



/* solve d^2/d(rho)^2 functions by direct iteration */
static xdouble iterdd_picard(sphr_t *sphr, xdouble rho,
    xdouble *ddcr, xdouble *ddtr, xdouble *ddck, xdouble *ddtk,
    const xdouble *dtr, const xdouble *dck, const xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    xdouble dmp, int itmax, xdouble tol)
{
  int it;
  xdouble err = 0;

  for ( it = 0; it < itmax; it++ ) {
    err = stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
        dtr, dck, dtk, tr, ck, tk, bphilr, fr, NULL, dmp);
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "iterdd_picard finished in %d steps, err %g\n", it, (double) err);
  return err;
}



/* advanced integral equation solvers */
const xdouble errinf = 1e300;
#include "ielmv.h"
#include "iemdiis.h"



/* load a previous solution for refinement */
static int loadcr(int npt, xdouble *ri, xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *ddcr, xdouble *ddtr,
    const char *fn)
{
  int i;
  FILE *fp;
  char ln[40960];
  xdouble r, c, t, dc, dt, ddc, ddt;

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

    if ( 7 != sscanf(ln,
          "%" XDBLSCNF "f %" XDBLSCNF "f %" XDBLSCNF "f "
          "%" XDBLSCNF "f %" XDBLSCNF "f"
          "%" XDBLSCNF "f %" XDBLSCNF "f",
          &r, &c, &t, &dc, &dt, &ddc, &ddt) ) {
      fprintf(stderr, "%s no data on line %d\n%s\n", fn, i, ln);
      break;
    }
    if ( FABS(r - ri[i]) > 1e-6 ) {
      fprintf(stderr, "%s: %d, %g %g mismatch\n",
          fn, i, (double) r, (double) ri[i]);
      break;
    }
    cr[i] = c;
    tr[i] = t;
    dcr[i] = dc;
    dtr[i] = dt;
    ddcr[i] = ddc;
    ddtr[i] = ddt;
  }
  fprintf(stderr, "%s loaded successfully\n", fn);
EXIT:
  fclose(fp);
  return i < npt;
}



static void savecr(int npt, xdouble *ri,
    xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr,
    xdouble *ddcr, xdouble *ddtr,
    xdouble *fr, xdouble *bphi, xdouble *bphilr, char *fn)
{
  int i;
  FILE *fp;

  xfopen(fp, fn, "w", return);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%8.6f %22.14e %22.14e %22.14e %22.14e %22.14e %22.14e %14.8f %g %g\n",
        (double) ri[i],
        (double) cr[i], (double) tr[i],
        (double) dcr[i], (double) dtr[i],
        (double) ddcr[i], (double) ddtr[i],
        (double) fr[i], (double) bphi[i], (double) bphilr[i]);
  }
  fclose(fp);
}



/* get pressure related variables */
static xdouble getpvars(sphr_t *sphr, xdouble rho,
    const xdouble *cr, const xdouble *tr,
    const xdouble *dcr, const xdouble *dtr,
    const xdouble *ddtr,
    const xdouble *bphilr, const xdouble *rdfr,
    xdouble *comprc, xdouble *dcomprc,
    xdouble *comprv, xdouble *dcomprv)
{
  int i;
  xdouble bpresv, yr, rDm1, dydt, ddydtt;

  *comprc = 1;
  *dcomprc = 0;
  bpresv = rho;
  *comprv = 1;
  *dcomprv = 0;
  for ( i = 0; i < sphr->npt; i++ ) {
    rDm1 = sphr->rDm1[i];
    /* comprc = d(beta P) / d(rho) = 1 - rho Int c(r) dr */
    *comprc -= rho * cr[i] * rDm1;
    /* dcomprc = d^2(beta P) / d(rho)^2 = -Int {c(r) + rho * d[c(r)]/ d(rho)} dr */
    *dcomprc -= (cr[i] + rho * dcr[i]) * rDm1;
    /* bpresv = rho + rho*rho/(2*D) Int r f'(r) y(r) dr */
    yr = getyr(tr[i], bphilr[i], ietype, params, &dydt, &ddydtt);
    bpresv += rho*rho/(2*dim) * rdfr[i] * yr * rDm1;
    /* d(bpresv)/d(rho) = 1 + rho/D Int r f'(r) y(r) dr
     *   + (rho*rho)/(2*D) Int r f'(r) d[y(r)]/d(rho) dr
     * where
     *  d[y(r)]/d(rho) = dydt d[t(r)]/d(rho) */
    *comprv += (rho/dim) * rdfr[i] * yr * rDm1
      + rho*rho/(2*dim) * rdfr[i] * dydt * dtr[i] * rDm1;
    /* d^2(bpresv)/d(rho)^2 = 1/D Int r f'(r) y(r) dr
     *   + 2*rho/D Int r f'(r) d[y(r)]/d(rho) dr
     *   + (rho*rho)/(2*D) Int r f'(r) d^2[y(r)]/d(rho)^2 dr
     * where
     *  d^2[y(r)]/d(rho)^2 = ddydtt { d[t(r)]/d(rho) }^2 + dydt d^2[t(r)]/d(rho)^2 */
    *dcomprv += (1./dim) * rdfr[i] * yr * rDm1
      + (2*rho/dim) * rdfr[i] * dydt * dtr[i] * rDm1
      + rho*rho/(2*dim) * rdfr[i] * (ddydtt * dtr[i] * dtr[i] + dydt * ddtr[i]) * rDm1;
  }
  return bpresv;
}



#if 0
/* compute the beta pressure */
static xdouble getpres_py(int npt, xdouble rho,
    xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *rdfr,
    xdouble *ck,
    xdouble *ri2, xdouble *ki2, xdouble *ri, xdouble *ki,
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
#endif



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
  xdouble *bphi, *bphisr, *bphilr, *fr, *rdfr, *fk;
  xdouble *cr, *ck, *tr, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
  xdouble *ddcr, *ddtr, *ddck, *ddtk;
  sphr_t *sphr;
  xdouble comprc, dcomprc;
  xdouble presv, comprv, dcomprv;
  xdouble B2;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);

  MAKE1DARR(bphi,   npt);
  MAKE1DARR(bphisr, npt);
  MAKE1DARR(bphilr, npt);
  MAKE1DARR(fr,     npt);
  MAKE1DARR(rdfr,   npt);
  MAKE1DARR(cr,     npt);
  MAKE1DARR(tr,     npt);
  MAKE1DARR(fk,     npt);
  MAKE1DARR(ck,     npt);
  MAKE1DARR(tk,     npt);
  MAKE1DARR(dcr,    npt);
  MAKE1DARR(dtr,    npt);
  MAKE1DARR(dck,    npt);
  MAKE1DARR(dtk,    npt);
  MAKE1DARR(ddcr,   npt);
  MAKE1DARR(ddtr,   npt);
  MAKE1DARR(ddck,   npt);
  MAKE1DARR(ddtk,   npt);

  mkfr(npt, beta, bphi, bphisr, bphilr,
      fr, rdfr, sphr->ri, sphr->dm, 0, 0, ljtype);
  sphr_r2k(sphr, fr, fk);
  COPY1DARR(cr, fr, npt);
  B2 = getB2(bphi, sphr->rDm1, sphr->npt);

  if ( fncrinp != NULL ) {
    if ( loadcr(npt, sphr->ri, cr, tr, dcr, dtr, ddcr, ddtr, fncrinp) != 0 )
      goto EXIT;
    savecr(npt, sphr->ri, cr, tr, dcr, dtr, ddcr, ddtr, fr, bphi, bphilr, "initcr.dat");
  }

  if ( solver == SOLVER_LMV ) {
    iter_lmv(sphr, rho, cr, tr, ck, tk,
        bphilr, fr, Mpt, damp, itmax, tol);
    iterd_lmv(sphr, rho, dcr, dtr, dck, dtk,
        tr, ck, tk, bphilr, fr, Mpt, damp, itmax, tol);
    iterdd_lmv(sphr, rho, ddcr, ddtr, ddck, ddtk,
        dtr, dck, dtk, tr, ck, tk, bphilr, fr, Mpt, damp, itmax, tol);
  } else if ( solver == SOLVER_MDIIS ) {
    iter_mdiis(sphr, rho, cr, tr, ck, tk,
        bphilr, fr, nbases, damp, itmax, tol);
    iterd_mdiis(sphr, rho, dcr, dtr, dck, dtk,
        tr, ck, tk, bphilr, fr, nbases, damp, itmax, tol);
    iterdd_mdiis(sphr, rho, ddcr, ddtr, ddck, ddtk,
        dtr, dck, dtk, tr, ck, tk, bphilr, fr, nbases, damp, itmax, tol);
  } else {
    iter_picard(sphr, rho, cr, tr, ck, tk,
        bphilr, fr, damp, itmax, tol);
    iterd_picard(sphr, rho, dcr, dtr, dck, dtk,
        tr, ck, tk, bphilr, fr, damp, itmax, tol);
    iterdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
        dtr, dck, dtk, tr, ck, tk, bphilr, fr, damp, itmax, tol);
  }

  if ( fncrout != NULL )
    savecr(npt, sphr->ri, cr, tr, dcr, dtr, ddcr, ddtr, fr, bphi, bphilr, fncrout);

  presv = getpvars(sphr, rho, cr, tr, dcr, dtr, ddtr, bphilr, rdfr,
      &comprc, &dcomprc, &comprv, &dcomprv);
  printf("rho %6.4f, T %7.4f: bPv %9.5f, dbPc/v %9.5f/%9.5f, ddbPc/v %9.5f/%9.5f, B2 %g\n",
      (double) rho, (double) (1/beta),
      (double) presv, (double) comprc, (double) comprv,
      (double) dcomprc, (double) dcomprv, (double) B2);
EXIT:
  sphr_close(sphr);
  FREE1DARR(bphi,   npt);
  FREE1DARR(bphisr, npt);
  FREE1DARR(bphilr, npt);
  FREE1DARR(fr,     npt);
  FREE1DARR(rdfr,   npt);
  FREE1DARR(cr,     npt);
  FREE1DARR(tr,     npt);
  FREE1DARR(fk,     npt);
  FREE1DARR(ck,     npt);
  FREE1DARR(tk,     npt);
  FREE1DARR(dcr,    npt);
  FREE1DARR(dtr,    npt);
  FREE1DARR(dck,    npt);
  FREE1DARR(dtk,    npt);
  FREE1DARR(ddcr,   npt);
  FREE1DARR(ddtr,   npt);
  FREE1DARR(ddck,   npt);
  FREE1DARR(ddtk,   npt);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  integ(numpt, rmax);
  return 0;
}

