#ifndef IELMV_H__
#define IELMV_H__



/* LMV solver
 * Reference:
 * Stanislav Labik, Anatol Malijevsky, Petr Vonka
 * A rapidly convergent method of solving the OZ equation
 * Molecular Physics, 1985, Vol. 56, No. 3, 709-715 */



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



static void iterlmv(sphr_t *sphr, xdouble rho,
    xdouble *fr, xdouble *der,
    xdouble *cr, xdouble *ck,
    xdouble *tr, xdouble *tk,
    xdouble *tr1, xdouble *tk1,
    int itmax, xdouble tol, xdouble dmp, int Mpt)
{
  int i, j, it, npt = sphr->npt, M;
  xdouble *Cjk = NULL, *mat = NULL, *a = NULL;
  xdouble *dp = NULL, *costab = NULL;
  xdouble c, dc, del, err1, err2, err, errp = 1e9;
  xdouble *crbest, errmin = 1e9;

  MAKE1DARR(crbest, npt);
  COPY1DARR(crbest, cr, npt);

  /* 1. initialize t(k) and t(r) */
  iterpicard(sphr, rho, cr, tr, ck, tk, fr, 0, 1);
  COPY1DARR(tr1, tr, npt);

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
    err = 0;
    for ( i = 0; i < npt; i++ ) {
      c = getcr(tr1[i], fr[i], ietype, params, &dc, NULL, NULL);
      del = c - cr[i];
      if ( FABS(del) > err ) err = FABS(del);
    }
    if ( err < errmin ) {
      COPY1DARR(crbest, cr, npt);
      errmin = err;
    }

    /* 2. compute c(r) and c(k) */
    for ( i = 0; i < npt; i++ ) {
      dc = 0;
      cr[i] = getcr(tr[i], fr[i], ietype, params, &dc, NULL, NULL);
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
      tk[j] += dmp * del;
    }

    /* 8. Use the OZ relation to solve for t(k) of large k */
    err2 = 0;
    for ( i = M; i < npt; i++ ) {
      del = rho * ck[i] * ck[i] / (1 - rho * ck[i]) - tk[i];
      if ( FABS(del) > err2 ) err2 = FABS(del);
      tk[i] += dmp * del;
    }
    sphr_k2r(sphr, tk, tr);

    for ( i = 0; i < npt; i++ )
      tk1[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    sphr_k2r(sphr, tk1, tr1);

    if (verbose)
      fprintf(stderr, "it %d: M %d, err %g -> %g, tkerr %g/%g, damp %g\n",
          it, M, (double) errp, (double) err,
          (double) err1, (double) err2, (double) dmp);
    if ( err < tol || err > errmax ) break;
    errp = err;
  }

  if ( verbose || err > tol ) {
    fprintf(stderr, "LMV stops after %d iterations, err %g/%g\n", it, err, errmin);
    if ( verbose >= 3 ) getchar();
  }

  /* use the best cr discovered so far */
  COPY1DARR(cr, crbest, npt);
  /* compute ck, tr, tk */
  iterpicard(sphr, rho, cr, tr, ck, tk, fr, 0, 1);

  FREE1DARR(crbest, npt);
  free(Cjk);
  free(mat);
  free(a);
  free(dp);
  free(costab);
}



#endif /* defined(LMV_H__) */

