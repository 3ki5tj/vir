/* integral equation */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



enum { SOLVER_PICARD, SOLVER_LMV, SOLVER_COUNT };



int dim = D;
int numpt = 8192;
int ffttype = 1;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 1.35;
xdouble beta;
xdouble rho = (xdouble) 0.278L;
int ietype = IETYPE_PY;
int solver = SOLVER_LMV;
int itmax = 10000;
xdouble tol = (xdouble) 1e-10L;
xdouble damp = 1;
const xdouble errmax = 1000;
int Mpt = 20;
char *fncrtr = NULL;
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
  argopt_add(ao, "-s", "%" XDBLSCNF "s", &params, "parameter of the integral equation");
  argopt_add(ao, "--damp", "%" XDBLSCNF "f", &damp, "damping factor of solving integral equations");
  argopt_add(ao, "--tol", "%" XDBLSCNF "f", &tol, "error tolerance");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file for cr and tr");
  argopt_add(ao, "-M", "%d", &Mpt, "number of points for LMV");
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



/* compute the cavity distribution function
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






/* direct (Picard iteration) */
static double iter(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *fr, int itmax)
{
  int i, it, npt = sphr->npt;
  xdouble x, dx, err;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, cr, ck);
    for ( i = 0; i < npt; i++ )
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    sphr_k2r(sphr, tk, tr);
    for ( err = 0, i = 0; i < npt; i++ ) {
      x = getcr(tr[i], fr[i], ietype, 0, NULL, NULL, NULL);
      dx = x - cr[i];
      if ( FABS(dx) > err ) err = FABS(dx);
      cr[i] += damp * dx;
    }
    if ( err < tol || err > errmax ) break;
    if ( verbose >= 2 ) fprintf(stderr, "it %d, err %g\n", it, (double) err);
  }
  if ( verbose || err > tol )
    fprintf(stderr, "main Picard finished in %d steps, err %g\n", it, (double) err);
  return err;
}



/* solve d/d(rho) functions by direct iteration */
static void iterd(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *ck, const xdouble *tk,
    xdouble *fr, int itmax)
{
  int i, it, npt = sphr->npt;
  xdouble x, dx, hk, err;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, dcr, dck);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dtk[i] = hk*hk + rho*hk*(2 + rho*hk)*dck[i];
    }
    sphr_k2r(sphr, dtk, dtr);
    for ( err = 0, i = 0; i < npt; i++ ) {
      x = fr[i] * dtr[i];
      dx = x - dcr[i];
      if ( FABS(dx) > err ) err = FABS(dx);
      dcr[i] += damp * dx;
    }
    if ( err < tol ) break;
  }
}



/* solve A x = b by L U decomposition
 * the matrix `a' will be destroyed
 * the vector `b' will be replaced by `x' on return */
static int lusolve(xdouble *a, xdouble *b, int n, xdouble tiny)
{
  int i, j, k, ip = 0;
  xdouble x, max;

  for (i = 0; i < n; i++) {  /* normalize each equation */
    for (max = 0.0, j = 0; j < n; j++)
      if ((x = FABS(a[i*n + j])) > max)
        max = x;
    if ( max <= 0 ) return -1;
    for (x = 1./max, j = 0; j < n; j++)
      a[i*n + j] *= x;
    b[i] *= x;
  }

  /* step 1: A = L U, column by column */
  for (j = 0; j < n; j++) {
    /* matrix U */
    for (i = 0; i < j; i++) {
      for (x = a[i*n + j], k = 0; k < i; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
    }

    /* matrix L, diagonal of L are 1 */
    max = 0.0;
    ip = j;
    for (i = j; i < n; i++) {
      for (x = a[i*n + j], k = 0; k < j; k++)
        x -= a[i*n + k] * a[k*n + j];
      a[i*n + j] = x;
      if ( FABS(x) > max ) {
        max = FABS(x);
        ip = i;
      }
    }

    if (j != ip) { /* swap the pivot row with the jth row */
      for (k = 0; k < n; k++)
        x = a[ip*n + k], a[ip*n + k] = a[j*n + k], a[j*n + k] = x;
      x = b[ip], b[ip] = b[j], b[j] = x;
    }
    if ( FABS(a[j*n + j]) < tiny )
      a[j*n + j] = tiny;
    /* divide by the pivot element, for the L matrix */
    if (j != n - 1)
      for (x = 1./a[j*n + j], i = j + 1; i < n; i++)
        a[i*n + j] *= x;
  }

  /* step 2: solve the equation L U x = b */
  for (i = 0; i < n; i++) { /* L y = b */
    x = b[i];
    for (j = 0; j < i; j++) x -= a[i*n + j] * b[j];
    b[i] = x;
  }
  for (i = n - 1; i >= 0; i--) { /* U x = y. */
    x = b[i];
    for (j = i + 1; j < n; j++) x -= a[i*n + j] * b[j];
    b[i] = x / a[i*n + i];
  }
  return 0;
}



/* compute Cjk */
static void getCjk(xdouble *Cjk, int npt, int M,
    xdouble *der, xdouble *costab, xdouble *dp)
{
  int m, k, l;

  for ( m = 1; m < 3*M - 1; m++ ) {
    for ( dp[m] = 0, l = 0; l < npt; l++ )
      dp[m] += der[l] * costab[m*npt + l];
    dp[m] /= npt;
  }

  for ( m = 0; m < M; m++ )
    for ( k = 0; k < M; k++ )
      Cjk[m*M + k] = dp[k - m + M] - dp[k + m + M];
}



/* compute the Jacobian matrix for the Newton-Raphson method */
static void getjacob(xdouble *mat, xdouble *b, int M,
    xdouble *ck, xdouble *tk, xdouble *Cjk,
    xdouble *ki, xdouble rho)
{
  int j, k;
  xdouble y;

  for ( j = 0; j < M; j++ ) {
    y = rho * ck[j] / (1 - rho * ck[j]);
    b[j] = ki[j] * (y * ck[j] - tk[j]);
    for ( k = 0; k < M; k++ )
      mat[j*M+k] = (k == j ? 1 : 0) - y * (2 + y) * Cjk[j*M+k];
  }
}



/* Reference:
 * Stanislav Labik, Anatol Malijevsky, Petr Vonka
 * A rapidly convergent method of solving the OZ equation
 * Molecular Physics, 1985, Vol. 56, No. 3, 709-715 */
static void iterlmv(sphr_t *sphr, xdouble rho,
    xdouble *fr, xdouble *fk, xdouble *der,
    xdouble *cr, xdouble *ck, xdouble *tr, xdouble *tk,
    int itmax, xdouble tol, int Mpt)
{
  int i, j, it, npt = sphr->npt, M;
  xdouble *Cjk = NULL, *mat = NULL, *a = NULL;
  xdouble *dp = NULL, *costab = NULL;
  xdouble c, dc, del, err1, err2, err, errp = 1e9;

  /* 1. initialize t(k) and t(r) */
  for ( i = 0; i < npt; i++ )
    tk[i] = rho * fk[i] * fk[i];
  sphr_k2r(sphr, tk, tr);

  /* 3. set M */
  if ( Mpt >= 0 ) {
    M = Mpt;
  } else { /* default value */
    M = (int) (4 * npt * sphr->dr);
  }
  if ( M >= npt ) M = npt;

  if ( verbose ) fprintf(stderr, "select M = %d\n", M);
  if ( M > 0 ) {
    xnew(Cjk, M*M);
    xnew(mat, M*M);
    xnew(a, M);
    xnew(dp, 3*M);
    xnew(costab, 3*M*npt);
    for ( j = 0; j < 3*M; j++ )
      for ( i = 0; i < npt; i++ )
        costab[j*npt + i] = COS(PI*(i*2+1)*(j-M)/npt/2);
  }

  for ( it = 0; it < itmax; it++ ) {
    /* 2. compute c(r) and c(k) */
    err = 0;
    for ( i = 0; i < npt; i++ ) {
      dc = 0;
      c = getcr(tr[i], fr[i], ietype, params, &dc, NULL, NULL);
      del = c - cr[i];
      if ( FABS(del) > err ) err = FABS(del);
      cr[i] += del;
      der[i] = dc;
    }
    sphr_r2k(sphr, cr, ck);

    /* 4. compute Cjk */
    getCjk(Cjk, npt, M, der, costab, dp);

    /* 6. compute the matrix for the Newton-Raphson method */
    getjacob(mat, a, M, ck, tk, Cjk, sphr->ki, rho);

    if ( lusolve(mat, a, M, 1e-10) != 0 )
      break;

    /* 7. Use the Newton-Raphson method to solve for t(k) of small k */
    err1 = 0;
    for ( j = 0; j < M; j++ ) {
      del = a[j] / sphr->ki[j];
      if ( FABS(del) > err1 ) err1 = FABS(del);
      tk[j] += damp * del;
    }

    /* 8. Use the OZ relation to solve for t(k) of large k */
    err2 = 0;
    for ( i = M; i < npt; i++ ) {
      del = rho * ck[i] * ck[i] / (1 - rho * ck[i]) - tk[i];
      if ( FABS(del) > err2 ) err2 = FABS(del);
      tk[i] += damp * del;
    }
    sphr_k2r(sphr, tk, tr);

    if (verbose)
      fprintf(stderr, "it %d: M %d, err %g -> %g, tkerr %g/%g, damp %g\n",
          it, M, (double) errp, (double) err, (double) err1, (double) err2, (double) damp);
    if ( err < tol || err > errmax ) break;
    //if (err > errp) {
    //  damp *= 0.8;
    //} else { /* try to recover to damp = 1 */
    //  damp = damp * 0.9 + 0.1;
    //}
    errp = err;
  }

  if ( verbose || err > tol )
    fprintf(stderr, "LMV stops after %d iterations, err %g\n", it, err);

  free(Cjk);
  free(mat);
  free(a);
  free(dp);
  free(costab);
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
  xdouble *cr, *tr, *ck, *tk;
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

  if ( solver == SOLVER_LMV ) {
    iterlmv(sphr, rho, fr, fk, der, cr, ck, tr, tk, itmax, tol, Mpt);
  } else {
    iter(sphr, rho, cr, tr, ck, tk, fr, itmax);
  }
  iterd(sphr, rho, dcr, dtr, dck, dtk, ck, tk, fr, itmax);
  if ( fncrtr != NULL )
    output(npt, sphr->ri, cr, tr, dcr, dtr, fr, bphi, fncrtr);
  pres = getpres_py(npt, rho, cr, tr, dcr, dtr, rdfr,
      ck, sphr->rDm1, sphr->kDm1, sphr->ri, sphr->ki,
      &compr, &dcompr, &presb, &comprb, &dcomprb, &presa);

  printf("rho %5.3f, T %6.3f: Pc %9.5f/%9.5f, dPc %9.5f, ddPc %9.5f, "
      "Pv %9.5f, dPv %9.5f, ddPv %9.5f\n",
      (double) rho, (double) (1/beta),
      (double) pres, (double) presa, (double) compr, (double) dcompr,
      (double) presb, (double) comprb, (double) dcomprb);

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

