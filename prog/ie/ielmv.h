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
  xdouble *der;
  xdouble *dek;
  xdouble *crbest;
  xdouble errmin;
  xdouble err1; /* error of tk[i] for i < M */
  xdouble err2; /* error of tk[i] for i >= M */
} lmv_t;



/* open an lmv object */
static lmv_t *lmv_open(int npt, int M, xdouble *ki)
{
  int i, j;
  lmv_t *lmv;

  xnew(lmv, 1);
  lmv->npt = npt;
  if ( M >= npt ) M = npt;
  if ( verbose ) fprintf(stderr, "LMV: select M = %d\n", M);
  lmv->M = M;
  lmv->ki = ki;
  MAKE1DARR(lmv->tr1,     npt);
  MAKE1DARR(lmv->tk1,     npt);
  MAKE1DARR(lmv->der,     npt);
  MAKE1DARR(lmv->dek,     npt);
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
  FREE1DARR(lmv->der,     lmv->npt);
  FREE1DARR(lmv->dek,     lmv->npt);
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



/* compute the Jacobian matrix for the Newton-Raphson method */
static void lmv_getjacob(lmv_t *lmv, xdouble *tk)
{
  int m, k, l, npt = lmv->npt, M = lmv->M;

  for ( m = 1; m < 3*M - 1; m++ ) {
    for ( lmv->dp[m] = 0, l = 0; l < npt; l++ )
      lmv->dp[m] += lmv->der[l] * lmv->costab[m*npt + l];
    lmv->dp[m] /= npt;
  }

  /* compute Cjk, d ck / dtk */
  for ( m = 0; m < M; m++ )
    for ( k = 0; k < M; k++ )
      lmv->Cjk[m*M + k] = lmv->dp[k - m + M] - lmv->dp[k + m + M];

  /* compute the Jacobian matrix */
  for ( m = 0; m < M; m++ ) {
    lmv->a[m] = lmv->ki[m] * (lmv->tk1[m] - tk[m]);
    for ( k = 0; k < M; k++ )
      lmv->mat[m*M+k] = (k == m ? 1 : 0) - lmv->dek[m] * lmv->Cjk[m*M+k];
  }
}



/* update tk */
static int lmv_update(lmv_t *lmv, xdouble *tk, xdouble dmp)
{
  int i, npt = lmv->npt, M = lmv->M;
  xdouble del;

  /* compute the Jacobian matrix for the Newton-Raphson method */
  lmv_getjacob(lmv, tk);

  if ( lusolve(lmv->mat, lmv->a, lmv->M, 1e-10) != 0 )
    return -1;

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

  return 0;
}



/* compute correlation functions using the LMV solver */
static xdouble iter_lmv(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int Mpt, xdouble dmp, int itmax, xdouble tol)
{
  int it, npt = sphr->npt;
  xdouble err, errp;
  lmv_t *lmv;

  /* open an lmv object */
  lmv = lmv_open(npt, Mpt, sphr->ki);
  COPY1DARR(lmv->crbest, cr, npt);

  /* initialize t(k) and t(r) */
  lmv->errmin = errp = step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, NULL, 0.0);
  COPY1DARR(lmv->tk1, tk, npt);

  for ( it = 0; it < itmax; it++ ) {
    /* compute the error of the current c(r) and c(k) */
    sphr_k2r(sphr, lmv->tk1, lmv->tr1);
    err = closure(cr, lmv->tr1, bphilr, fr, NULL, NULL,
        npt, ietype, params, 0.0);
    lmv_savebest(lmv, cr, err);

    /* compute c(r) and c(k) */
    closure(cr, tr, bphilr, fr, NULL, lmv->der,
        npt, ietype, params, 1.0);
    sphr_r2k(sphr, cr, ck);
    oz(ck, lmv->tk1, lmv->dek, npt, rho);

    /* compute the new tk */
    if ( lmv_update(lmv, tk, dmp) != 0 ) break;
    sphr_k2r(sphr, tk, tr);

    if ( verbose >= 2 )
      fprintf(stderr, "iter_lmv: it %d: M %d, err %g -> %g (min %g), tkerr %g/%g, damp %g\n",
          it, lmv->M, (double) errp, (double) err, (double) lmv->errmin,
          (double) lmv->err1, (double) lmv->err2, (double) dmp);
    if ( err < tol || err > errmax ) break;
    errp = err;
  }

  /* use the best cr discovered so far */
  COPY1DARR(cr, lmv->crbest, npt);
  /* compute ck, tr, tk */
  err = step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, NULL, 0.0);
  if ( verbose )
    fprintf(stderr, "iter_lmv finishes in %d iterations, err %g\n", it, (double) err);

  lmv_close(lmv);
  return err;
}



/* compute d/d(rho) functions using the LMV solver */
static xdouble iterd_lmv(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int Mpt, xdouble dmp, int itmax, xdouble tol)
{
  int it, npt = sphr->npt;
  xdouble err, errp;
  lmv_t *lmv;

  /* open an lmv object */
  lmv = lmv_open(npt, Mpt, sphr->ki);
  COPY1DARR(lmv->crbest, dcr, npt);

  /* initialize t(k) and t(r) */
  lmv->errmin = errp = stepd_picard(sphr, rho, dcr, dtr, dck, dtk,
      tr, ck, tk, bphilr, fr, NULL, 0.0);
  COPY1DARR(lmv->tk1, dtk, npt);

  for ( it = 0; it < itmax; it++ ) {
    /* compute the error of the current c(r) and c(k) */
    sphr_k2r(sphr, lmv->tk1, lmv->tr1);
    err = closured(dcr, lmv->tr1, tr, bphilr, fr, NULL, NULL,
        npt, ietype, params, 0.0);
    lmv_savebest(lmv, dcr, err);

    /* compute dc(r) and dc(k) */
    closured(dcr, dtr, tr, bphilr, fr, NULL, lmv->der,
        npt, ietype, params, 1.0);
    sphr_r2k(sphr, dcr, dck);
    ozd(dck, lmv->tk1, ck, tk, lmv->dek, npt, rho);

    /* compute the new dtk */
    if ( lmv_update(lmv, dtk, dmp) != 0 ) break;
    sphr_k2r(sphr, dtk, dtr);

    if ( verbose >= 2 )
      fprintf(stderr, "iterd_lmv: it %d: M %d, err %g -> %g (min %g), tkerr %g/%g, damp %g\n",
          it, lmv->M, (double) errp, (double) err, (double) lmv->errmin,
          (double) lmv->err1, (double) lmv->err2, (double) dmp);
    if ( err < tol || err > errmax ) break;
    errp = err;
  }

  /* use the best dcr discovered so far */
  COPY1DARR(dcr, lmv->crbest, npt);
  /* compute dck, dtr, dtk */
  err = stepd_picard(sphr, rho, dcr, dtr, dck, dtk,
      tr, ck, tk, bphilr, fr, NULL, 0.0);
  if ( verbose || err > tol )
    fprintf(stderr, "iterd_lmv finishes in %d iterations, err %g\n", it, (double) err);

  lmv_close(lmv);
  return err;
}



/* compute d^2/d(rho)^2 functions using the LMV solver */
static xdouble iterdd_lmv(sphr_t *sphr, xdouble rho,
    xdouble *ddcr, xdouble *ddtr, xdouble *ddck, xdouble *ddtk,
    const xdouble *dtr, const xdouble *dck, const xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int Mpt, xdouble dmp, int itmax, xdouble tol)
{
  int it, npt = sphr->npt;
  xdouble err, errp;
  lmv_t *lmv;

  /* open an lmv object */
  lmv = lmv_open(npt, Mpt, sphr->ki);
  COPY1DARR(lmv->crbest, ddcr, npt);

  /* initialize t(k) and t(r) */
  lmv->errmin = errp = stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
      dtr, dck, dtk, tr, ck, tk, bphilr, fr, NULL, 0.0);
  COPY1DARR(lmv->tk1, ddtk, npt);

  for ( it = 0; it < itmax; it++ ) {
    /* compute the error of the current ddc(r) and ddc(k) */
    sphr_k2r(sphr, lmv->tk1, lmv->tr1);
    err = closuredd(ddcr, lmv->tr1, dtr, tr, bphilr, fr, NULL, NULL,
        npt, ietype, params, 0.0);
    lmv_savebest(lmv, ddcr, err);

    /* compute ddc(r) and ddc(k) */
    closuredd(ddcr, ddtr, dtr, tr, bphilr, fr, NULL, lmv->der,
        npt, ietype, params, 1.0);
    sphr_r2k(sphr, ddcr, ddck);
    ozdd(ddck, lmv->tk1, dck, dtk, ck, tk, lmv->dek, npt, rho);

    /* compute the new dtk */
    if ( lmv_update(lmv, ddtk, dmp) != 0 ) break;
    sphr_k2r(sphr, ddtk, ddtr);

    if ( verbose >= 2 )
      fprintf(stderr, "iterd_lmv: it %d: M %d, err %g -> %g (min %g), tkerr %g/%g, damp %g\n",
          it, lmv->M, (double) errp, (double) err, (double) lmv->errmin,
          (double) lmv->err1, (double) lmv->err2, (double) dmp);
    if ( err < tol || err > errmax ) break;
    errp = err;
  }

  /* use the best dcr discovered so far */
  COPY1DARR(ddcr, lmv->crbest, npt);
  /* compute dck, dtr, dtk */
  err = stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
      dtr, dck, dtk, tr, ck, tk, bphilr, fr, NULL, 0.0);
  if ( verbose || err > tol )
    fprintf(stderr, "iterd_lmv finishes in %d iterations, err %g\n", it, (double) err);

  lmv_close(lmv);
  return err;
}


#endif /* defined(LMV_H__) */

