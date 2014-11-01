#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"


int dim = D;
int numpt = 1024;
int ffttype = 1;
xdouble rmax = (xdouble) 20.48L;
xdouble rhomax = (xdouble) 0.8L;
int itmax = 10000;
xdouble tol = (xdouble) 1e-8L;
int Mpt = -1;
char *fncrtr = NULL;
int verbose = 0;


/*
enum {
  IETYPE_PY = 0,
  IETYPE_HNC = 1,
  IETYPE_IR = 2,
  IETYPE_HC = 3,
  IETYPE_SQR = 4
};
*/


int ietype = IETYPE_PY;
xdouble ie_s = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  int dohnc = 0;

  ao->desc = "solve integral equation";
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rhomax, "maximal rho");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--itmax", "%d", &itmax, "maximal number of iterations");
  argopt_add(ao, "--tol", "%" XDBLSCNF "f", &tol, "tolerance");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file for saving cr and tr");
  argopt_add(ao, "--Mpt", "%d", &Mpt, "number of points for the Newton-Raphson method");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted-chain approximation");
  argopt_add(ao, "-s", "%" XDBLSCNF "f", &ie_s, "s of the integral equation");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( dohnc ) ietype = IETYPE_HNC;
  printf("rmax %f, rho %g, ietype %d\n",
      (double) rmax, (double) rhomax, ietype);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute the cavity function
 * *dy = partial y / partial t  */
__inline static xdouble getyr(xdouble tr, int ietype, xdouble s,
    xdouble *dy, xdouble *w, xdouble *dw)
{
  xdouble u, z;

  if ( ietype == IETYPE_PY ) { /* PY */
    if (dy) *dy = 1;
    return 1 + tr;
  } else if ( ietype == IETYPE_HNC ) { /* HNC */
    z = EXP(tr);
    if (dy) *dy = z;
    if (w) *w = z - 1 - tr;
    if (dw) *dw = z - 1;
    return z;
  } else if ( ietype == IETYPE_INVROWLINSON ) { /* inverse Rowlinson */
    z = EXP(tr);
    if (dy) *dy = 1 + s * (z - 1);
    if (w) *w = z - 1 - tr;
    if (dw) *dw = z - 1;
    return 1 + tr + s * (z - 1 - tr);
  } else if ( ietype == IETYPE_HC ) { /* Hutchinson-Conkie */
    u = 1 + s * tr;
    if (u < XDBL_MIN) u = XDBL_MIN;
    z = POW(u, 1./s);
    if (dy) *dy = z/u;
    if (w) *w = (-LOG(u) + 1 - 1/u) * z / (s * s);
    if (dw) *dw = (-LOG(u) + (1 - s) * (1 - 1/u)) * z / u / (s * s);
    return z;
  } else { /* quadratic */
    if (dy) *dy = 1 + s * tr;
    if (w) *w = .5 * tr * tr;
    if (dw) *dw = tr;
    return 1 + tr + s * .5 * tr * tr;
  }
}



/* solve the linear equation: m x = b */
static int linsolve(int n, xdouble *m, xdouble *x, xdouble *b)
{
  int i, j, k, jm;
  xdouble y;

  for ( i = 0; i < n; i++ ) x[i] = 0;

  for ( i = 0; i < n; i++ ) {
    /* 1. select the pivot of the ith column
     * pivot: the maximal element m(j = i..n-1, i)  */
    y = FABS(m[(jm = i)*n + i]);
    for ( j = i + 1; j < n; j++ )
      if ( FABS(m[j*n + i]) > y )
        y = FABS(m[(jm = j)*n + i]);

    /* 2. swap the jm'th (pivot) row with the ith row */
    if ( jm != i ) {
      XDBL_SWAP(b[jm], b[i]);
      for ( j = i; j < n; j++ )
        XDBL_SWAP(m[i*n+j], m[jm*n+j]);
    }
    y = m[i*n+i];
    if ( FABS(y) < XDBL_MIN ) {
      fprintf(stderr, "singular matrix on %dth row\n", i);
      return -1;
    }

    /* 3. normalize the ith row */
    b[i] /= y;
    for ( k = i; k < n; k++ )
      m[i*n+k] /= y;

    /* 4. use the pivot row to eliminate the following rows */
    for ( j = i + 1; j < n; j++ ) { /* for rows */
      y = m[j*n+i];
      b[j] -= y * b[i];
      for ( k = i; k < n; k++ )
        m[j*n+k] -= y * m[i*n+k];
    }
  }

  /* 5. now that the matrix is upper-triangular
   *    solve for x */
  for ( i = n - 1; i >= 0; i-- ) {
    x[i] = b[i] / m[i*n+i];
    for ( j = 0; j < i; j++ ) {
      b[j] -= m[j*n+i] * x[i];
    }
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
    int ietype, xdouble s,
    xdouble *fr, xdouble *fk, xdouble *der,
    xdouble *cr, xdouble *ck, xdouble *tr, xdouble *tk,
    int itmax, xdouble tol, int Mpt)
{
  int i, j, it, npt = sphr->npt, M;
  xdouble *Cjk = NULL, *mat = NULL, *a = NULL, *b = NULL;
  xdouble *dp = NULL, *costab = NULL;
  xdouble y, dy, damp = 1, err1, err2, err, errp = 1e9;

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
    xnew(b, M);
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
      dy = 0;
      y = getyr(tr[i], ietype, s, &dy, NULL, NULL);
      y = (fr[i] + 1) * y - cr[i] - tr[i] - 1;
      if ( FABS(y) > err ) err = FABS(y);
      cr[i] += y;
      der[i] = (fr[i] + 1) * dy - 1;
    }
    sphr_r2k(sphr, cr, ck);

    /* 4. compute Cjk */
    getCjk(Cjk, npt, M, der, costab, dp);

    /* 6. compute the matrix for the Newton-Raphson method */
    getjacob(mat, b, M, ck, tk, Cjk, sphr->ki, rho);

    if ( linsolve(M, mat, a, b) != 0 )
      break;

    /* 7. Use the Newton-Raphson method to solve for t(k) of small k */
    err1 = 0;
    for ( j = 0; j < M; j++ ) {
      y = a[j] / sphr->ki[j];
      if ( FABS(y) > err1 ) err1 = FABS(y);
      tk[j] += damp * y;
    }

    /* 8. Use the OZ relation to solve for t(k) of large k */
    err2 = 0;
    for ( i = M; i < npt; i++ ) {
      y = rho * ck[i] * ck[i] / (1 - rho * ck[i]) - tk[i];
      if ( FABS(y) > err2 ) err2 = FABS(y);
      tk[i] += damp * y;
    }
    sphr_k2r(sphr, tk, tr);

    fprintf(stderr, "it %d: M %d, err %g -> %g, tkerr %g/%g, damp %g\n",
        it, M, (double) errp, (double) err, (double) err1, (double) err2, (double) damp);
    if (err < tol) break;
    if (err > errp) {
      damp *= 0.8;
    } else { /* try to recover to damp = 1 */
      damp = damp * 0.9 + 0.1;
    }
    errp = err;
  }

  if ( verbose || it >= itmax )
    fprintf(stderr, "iteration stops after %d iterations\n", it);

  free(Cjk);
  free(mat);
  free(a);
  free(b);
  free(dp);
  free(costab);
}



static int savecrtr1(const char *fn, sphr_t *sphr,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk)
{
  FILE *fp;
  int i, npt = sphr->npt;

  xfopen(fp, fn, "w", return -1);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%g %g %g %g %g %g\n",
        (double) sphr->ri[i], (double) cr[i], (double) tr[i],
        (double) sphr->ki[i], (double) ck[i], (double) tk[i]);
  }
  fclose(fp);
  return 0;
}



static void integ(int npt, xdouble rmax, xdouble rhomax)
{
  xdouble rho = rhomax;
  xdouble *bphi, *fr, *fk, *der, *cr, *ck, *tr, *tk;
  int i;
  sphr_t *sphr;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);
  xnew(bphi, npt);
  xnew(fr, npt);
  xnew(fk, npt);
  xnew(der, npt);
  xnew(cr, npt);
  xnew(ck, npt);
  xnew(tr, npt);
  xnew(tk, npt);

  /* solvent-solvent interaction */
  for ( i = 0; i < npt; i++ )
    fr[i] = (sphr->ri[i] < 1) ? -1 : 0;
  sphr_r2k(sphr, fr, fk);

  iterlmv(sphr, rho, ietype, ie_s, fr, fk, der,
      cr, ck, tr, tk, itmax, tol, Mpt);
  if ( fncrtr != NULL )
    savecrtr1(fncrtr, sphr, cr, tr, ck, tk);

  sphr_close(sphr);
  free(bphi);
  free(fr);
  free(fk);
  free(der);
  free(cr);
  free(ck);
  free(tr);
  free(tk);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  integ(numpt, rmax, rhomax);
  return 0;
}

