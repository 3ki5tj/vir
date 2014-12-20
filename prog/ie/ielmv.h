#ifndef IELMV_H__
#define IELMV_H__



/* LMV solver
 * Reference:
 * Stanislav Labik, Anatol Malijevsky, Petr Vonka
 * A rapidly convergent method of solving the OZ equation
 * Molecular Physics, 1985, Vol. 56, No. 3, 709-715 */



typedef struct {
  int npt;
  int M;
  xdouble *ki;
  xdouble *Cjk; /* d ck / d tk */
  xdouble *mat; /* M x M matrix */
  xdouble *a; /* M array */
  xdouble *dp; /* 3*M array */
  xdouble *costab; /* 3*M*npt */
  xdouble *tr1;
  xdouble *tk1;
  xdouble *crbest;
  xdouble errmin;
  xdouble err1; /* error of tk[i] for i < M */
  xdouble err2; /* error of tk[i] for i >= M */
} lmv_t;



/* open an lmv object */
static lmv_t *lmv_open(int npt, int M, xdouble dr, xdouble *ki)
{
  int i, j;
  lmv_t *lmv;

  if ( M < 0 ) M = (int) (4 * npt * dr);
  if ( M >= npt ) M = npt;
  if ( verbose ) fprintf(stderr, "select M = %d\n", M);

  xnew(lmv, 1);
  lmv->npt = npt;
  lmv->M = M;
  lmv->ki = ki;
  MAKE1DARR(lmv->tr1,     npt);
  MAKE1DARR(lmv->tk1,     npt);
  MAKE1DARR(lmv->crbest,  npt);
  lmv->errmin = errinf;

  if ( M > 0 ) {
    xnew(lmv->Cjk, M*M);
    xnew(lmv->mat, M*M);
    xnew(lmv->a,  M);
    xnew(lmv->dp, 3*M);
    xnew(lmv->costab, 3*M*npt);
    for ( j = 0; j < 3*M; j++ )
      for ( i = 0; i < npt; i++ )
        lmv->costab[j*npt + i] = COS(PI*(i*2+1)*(j-M)/npt/2);
  }

  return lmv;
}



static void lmv_close(lmv_t *lmv)
{
  if ( lmv == NULL ) return;
  FREE1DARR(lmv->tr1,     lmv->npt);
  FREE1DARR(lmv->tk1,     lmv->npt);
  FREE1DARR(lmv->crbest,  lmv->npt);
  if ( lmv->M > 0 ) {
    free(lmv->Cjk);
    free(lmv->mat);
    free(lmv->a);
    free(lmv->dp);
    free(lmv->costab);
  }
  free(lmv);
}



/* register a good cr */
static void lmv_savebest(lmv_t *lmv, xdouble *cr, xdouble err)
{
  if ( err < lmv->errmin ) {
    COPY1DARR(lmv->crbest, cr, lmv->npt);
    lmv->errmin = err;
  }
}



/* compute Cjk */
static void lmv_getCjk(lmv_t *lmv, xdouble *der)
{
  int m, k, l, npt = lmv->npt, M = lmv->M;

  for ( m = 1; m < 3*M - 1; m++ ) {
    for ( lmv->dp[m] = 0, l = 0; l < npt; l++ )
      lmv->dp[m] += der[l] * lmv->costab[m*npt + l];
    lmv->dp[m] /= npt;
  }

  for ( m = 0; m < M; m++ )
    for ( k = 0; k < M; k++ )
      lmv->Cjk[m*M + k] = lmv->dp[k - m + M] - lmv->dp[k + m + M];
}



/* compute the Jacobian matrix for the Newton-Raphson method */
static void lmv_getjacob(lmv_t *lmv, xdouble *ck, xdouble *tk,
    xdouble rho)
{
  int j, k, M = lmv->M;
  xdouble y;

  for ( j = 0; j < M; j++ ) {
    y = rho * ck[j] / (1 - rho * ck[j]);
    lmv->a[j] = lmv->ki[j] * (y * ck[j] - tk[j]);
    for ( k = 0; k < M; k++ )
      lmv->mat[j*M+k] = (k == j ? 1 : 0) - y * (2 + y) * lmv->Cjk[j*M+k];
  }
}



/* update tk */
static void lmv_update(lmv_t *lmv, xdouble *tk, xdouble dmp)
{
  int i, npt = lmv->npt, M = lmv->M;
  xdouble del;

  /* use the Newton-Raphson method to solve for t(k) of small k */
  lmv->err1 = 0;
  for ( i = 0; i < M; i++ ) {
    del = lmv->a[i] / lmv->ki[i];
    if ( FABS(del) > lmv->err1 ) lmv->err1 = FABS(del);
    tk[i] += dmp * del;
  }

  /* use the OZ relation to solve for t(k) of large k */
  lmv->err2 = 0;
  for ( i = M; i < npt; i++ ) {
    del = lmv->tk1[i] - tk[i];
    if ( FABS(del) > lmv->err2 ) lmv->err2 = FABS(del);
    tk[i] += dmp * del;
  }
}



static void iter_lmv(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *ck,
    xdouble *tr, xdouble *tk,
    xdouble *fr, xdouble *der,
    int Mpt, xdouble dmp, int itmax, xdouble tol)
{
  int it, npt = sphr->npt;
  xdouble err, errp = errinf;
  lmv_t *lmv;

  /* open an lmv object */
  lmv = lmv_open(npt, Mpt, sphr->dr, sphr->ki);
  COPY1DARR(lmv->crbest, cr, npt);

  /* initialize t(k) and t(r) */
  step_picard(sphr, rho, cr, tr, ck, tk, fr, NULL, 0.0);
  COPY1DARR(lmv->tk1, tk, npt);

  for ( it = 0; it < itmax; it++ ) {
    /* compute the error of the current c(r) and c(k) */
    sphr_k2r(sphr, lmv->tk1, lmv->tr1);
    err = closure(cr, lmv->tr1, fr, NULL, NULL, npt, ietype, params, 0.0);
    lmv_savebest(lmv, cr, err);

    /* compute c(r) and c(k) */
    closure(cr, tr, fr, NULL, der, npt, ietype, params, 1.0);
    sphr_r2k(sphr, cr, ck);

    /* compute Cjk */
    lmv_getCjk(lmv, der);

    /* compute the matrix for the Newton-Raphson method */
    lmv_getjacob(lmv, ck, tk, rho);

    if ( lusolve(lmv->mat, lmv->a, lmv->M, 1e-10) != 0 )
      break;

    oz(ck, lmv->tk1, npt, rho);
    lmv_update(lmv, tk, dmp);
    sphr_k2r(sphr, tk, tr);

    if (verbose)
      fprintf(stderr, "it %d: M %d, err %g -> %g, tkerr %g/%g, damp %g\n",
          it, lmv->M, (double) errp, (double) err,
          (double) lmv->err1, (double) lmv->err2, (double) dmp);
    if ( err < tol || err > errmax ) break;
    errp = err;
  }

  if ( verbose || err > tol ) {
    fprintf(stderr, "LMV stops after %d iterations, err %g/%g\n", it, err, lmv->errmin);
    if ( verbose >= 3 ) getchar();
  }

  /* use the best cr discovered so far */
  COPY1DARR(cr, lmv->crbest, npt);
  /* compute ck, tr, tk */
  step_picard(sphr, rho, cr, tr, ck, tk, fr, NULL, 0.0);
  lmv_close(lmv);
}



#endif /* defined(LMV_H__) */

