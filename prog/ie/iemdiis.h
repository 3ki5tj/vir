#ifndef IEMDIIS_H__
#define IEMDIIS_H__



/* MDIIS solver */



typedef struct {
  int npt;
  int mnb; /* maximal number of bases */
  int nb; /* number of functions in the basis */
  xdouble **cr; /* basis */
  xdouble **res; /* residues */
  xdouble *mat; /* correlations of residues */
  xdouble *mat2; /* temporary matrix for the LU decomposition */
  xdouble *coef; /* coefficients */
  xdouble *crbest;
  xdouble errmin;
} mdiis_t;



/* open an mdiis object */
static mdiis_t *mdiis_open(int npt, int mnb)
{
  mdiis_t *m;
  int mnb1;

  xnew(m, 1);
  m->npt = npt;
  m->mnb = mnb;
  m->nb = 0;
  mnb1 = mnb + 1;
  MAKE2DARR(m->cr, mnb1, npt);
  MAKE2DARR(m->res, mnb1, npt);
  MAKE1DARR(m->mat, mnb1 * mnb1);
  MAKE1DARR(m->mat2, mnb1 * mnb1);
  MAKE1DARR(m->coef, mnb1);
  MAKE1DARR(m->crbest, npt);
  m->errmin = errinf;
  return m;
}



/* close the mdiis object */
static void mdiis_close(mdiis_t *m)
{
  int mnb1, npt;
  if ( m == NULL ) return;
  mnb1 = m->mnb + 1;
  npt = m->npt;
  FREE2DARR(m->cr,      mnb1, npt);
  FREE2DARR(m->res,     mnb1, npt);
  FREE1DARR(m->mat,     mnb1 * mnb1);
  FREE1DARR(m->mat2,    mnb1 * mnb1);
  FREE1DARR(m->coef,    mnb1);
  FREE1DARR(m->crbest,  npt);
  free(m);
}



static int mdiis_solve(mdiis_t *m)
{
  int nb = m->nb, nb1 = m->nb + 1, mnb1 = m->mnb + 1, i, j;

  for ( i = 0; i < nb; i++ ) m->coef[i] = 0;
  m->coef[nb] = -1;
  /* copy the matrix, for the content is to be destroyed */
  for ( i = 0; i < nb1; i++ )
    for ( j = 0; j < nb1; j++ )
      m->mat2[i*nb1 + j] = m->mat[i*mnb1 + j];
  for ( i = 0; i < nb1; i++ )
    m->mat2[i*nb1 + nb] = m->mat2[nb*nb1 + i] = -1;
  m->mat2[nb*nb1 + nb] = 0;
  if ( lusolve(m->mat2, m->coef, nb1, 1e-20) != 0 ) {
    fprintf(stderr, "MDIIS lusolve failed\n");
    exit(1);
  }
  return 0;
}



/* construct the new c(r) */
static void mdiis_gencr(mdiis_t *m, xdouble *cr, xdouble damp)
{
  int npt = m->npt, nb = m->nb;
  int k, il;

  for ( il = 0; il < npt; il++ )
    m->cr[nb][il] = 0;
  for ( k = 0; k < nb; k++ ) {
    xdouble coef = m->coef[k];
    for ( il = 0; il < npt; il++ )
      m->cr[nb][il] += coef * (m->cr[k][il] + damp * m->res[k][il]);
  }

  /* save m->cr[nb] to c(r) */
  COPY1DARR(cr, m->cr[nb], npt);
}



/* compute the dot product */
static xdouble mdiis_getdot(xdouble *a, xdouble *b, int n)
{
  int i;
  xdouble x = 0;

  for ( i = 0; i < n; i++ ) x += a[i] * b[i];
  return x / n;
}



/* build the residue correlation matrix */
static int mdiis_build(mdiis_t *m, xdouble *cr, xdouble *res)
{
  int ib, mnb, mnb1, npt = m->npt;

  m->nb = 1;
  mnb = m->mnb;
  mnb1 = m->mnb + 1;

  COPY1DARR(m->cr[0], cr, npt);
  COPY1DARR(m->res[0], res, npt);

  m->mat[0] = mdiis_getdot(m->res[0], m->res[0], npt);
  for ( ib = 0; ib < mnb; ib++ )
    m->mat[ib*mnb1 + mnb] = m->mat[mnb*mnb1 + ib] = -1;
  m->mat[mnb*mnb1 + mnb] = 0;
  return 0;
}



/* replace base ib by cr */
static int mdiis_update(mdiis_t *m, xdouble *cr, xdouble *res, xdouble err)
{
  int i, ib, nb, mnb1, npt = m->npt;
  xdouble dot, max;

  nb = m->nb;
  mnb1 = m->mnb + 1;

  /* save this function if it achieves the minimal error so far */
  if ( err < m->errmin ) {
    COPY1DARR(m->crbest, m->cr[nb], npt);
    m->errmin = err;
  }

  if ( nb < m->mnb ) {
    ib = nb;
    m->nb = ++nb;
  } else {
    /* choose the base with the largest residue */
    ib = 0;
    for ( i = 1; i < nb; i++ )
      /* the diagonal represents the error */
      if ( m->mat[i*mnb1+i] > m->mat[ib*mnb1 + ib] )
        ib = i;
    max = m->mat[ib*mnb1 + ib];

    dot = mdiis_getdot(res, res, npt);
    if ( dot > max ) {
#ifndef MDIIS_THRESHOLD
#define MDIIS_THRESHOLD 10.0
#endif
      int reset = ( SQRT(dot) < MDIIS_THRESHOLD );
      if ( verbose ) {
        fprintf(stderr, "MDIIS: bad basis, %g is greater than %g, %s, error:",
          (double) dot, (double) max, reset ? "reset" : "accept");
        for ( i = 0; i < nb; i++ )
          fprintf(stderr, " %g", (double) m->mat[i*mnb1+i]);
        fprintf(stderr, "\n");
      }
      if ( reset ) {
        mdiis_build(m, cr, res);
        return 1;
      }
    }
  }

  /* replace base ib by cr */
  COPY1DARR(m->cr[ib], cr, npt);
  COPY1DARR(m->res[ib], res, npt);

  /* update the residue correlation matrix
   * note: we do not need to update the last row & column */
  for ( i = 0; i < nb; i++ )
    m->mat[i*mnb1 + ib] = m->mat[ib*mnb1 + i]
      = mdiis_getdot(m->res[i], res, npt);
  return ib;
}



/* compute correlation functions by the MDIIS method */
static xdouble iter_mdiis(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int nbases, xdouble dmp, int itmax, xdouble tol)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, npt = sphr->npt;
  xdouble err, errp, *res;

  /* open an mdiis object */
  mdiis = mdiis_open(npt, nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial basis */
  COPY1DARR(mdiis->crbest, cr, npt);
  step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, res, 0);
  mdiis_build(mdiis, cr, res);

  err = errp = errinf;
  for ( it = 0; it < itmax; it++ ) {
    /* compute a set of coefficients of combining basic functions */
    mdiis_solve(mdiis);
    /* construct a new cr from the combination */
    mdiis_gencr(mdiis, cr, dmp);
    /* compute the residue vector and error */
    err = step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, res, 0);
    /* add the new cr into the basis */
    ib = mdiis_update(mdiis, cr, res, err);

    if ( verbose >= 2 )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n",
          it, (double) errp, (double) err, ibp, ib);
    if ( err < tol ) break;
    ibp = ib;
    errp = err;
  }
  /* use the best cr discovered so far */
  COPY1DARR(cr, mdiis->crbest, npt);
  /* update the corresponding ck, tr, tk */
  err = step_picard(sphr, rho, cr, tr, ck, tk, bphilr, fr, NULL, 0);
  if ( verbose || err > tol )
    fprintf(stderr, "iter_mdiis finished in %d steps, err %g\n", it, (double) err);

  mdiis_close(mdiis);
  return err;
}



/* compute d/d(rho) functions by the MDIIS method */
static xdouble iterd_mdiis(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int nbases, xdouble dmp, int itmax, xdouble tol)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, npt = sphr->npt;
  xdouble err, errp, *res;

  /* open an mdiis object */
  mdiis = mdiis_open(npt, nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial basis */
  COPY1DARR(mdiis->crbest, dcr, npt);
  stepd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, bphilr, fr, res, 0);
  mdiis_build(mdiis, dcr, res);

  err = errp = errinf;
  for ( it = 0; it < itmax; it++ ) {
    /* compute a set of coefficients of combining basic functions */
    mdiis_solve(mdiis);
    /* construct a new dcr from the combination */
    mdiis_gencr(mdiis, dcr, dmp);
    /* compute the residue vector and error */
    err = stepd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, bphilr, fr, res, 0);
    /* add the new dcr into the basis */
    ib = mdiis_update(mdiis, dcr, res, err);

    if ( verbose >= 2 )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n",
          it, (double) errp, (double) err, ibp, ib);
    if ( err < tol ) break;
    ibp = ib;
    errp = err;
  }
  /* use the best cr discovered so far */
  COPY1DARR(dcr, mdiis->crbest, npt);
  /* update the corresponding ck, tr, tk */
  err = stepd_picard(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, bphilr, fr, NULL, 0);
  if ( verbose || err > tol )
    fprintf(stderr, "iterd_mdiis finished in %d steps, err %g\n", it, (double) err);

  mdiis_close(mdiis);
  return err;
}



/* compute d^2/d(rho)^2 functions by the MDIIS method */
static xdouble iterdd_mdiis(sphr_t *sphr, xdouble rho,
    xdouble *ddcr, xdouble *ddtr, xdouble *ddck, xdouble *ddtk,
    const xdouble *dtr, const xdouble *dck, const xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *bphilr, const xdouble *fr,
    int nbases, xdouble dmp, int itmax, xdouble tol)
{
  mdiis_t *mdiis;
  int it, ibp = 0, ib, npt = sphr->npt;
  xdouble err, errp, *res;

  /* open an mdiis object */
  mdiis = mdiis_open(npt, nbases);
  /* use the space of the last array for `res' */
  res = mdiis->res[mdiis->mnb];

  /* construct the initial basis */
  COPY1DARR(mdiis->crbest, ddcr, npt);
  stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
      dtr, dck, dtk, tr, ck, tk, bphilr, fr, res, 0);
  mdiis_build(mdiis, ddcr, res);

  err = errp = errinf;
  for ( it = 0; it < itmax; it++ ) {
    /* compute a set of coefficients of combining basic functions */
    mdiis_solve(mdiis);
    /* construct a new ddcr from the combination */
    mdiis_gencr(mdiis, ddcr, dmp);
    /* compute the residue vector and error */
    err = stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
        dtr, dck, dtk, tr, ck, tk, bphilr, fr, res, 0);
    /* add the new ddcr into the basis */
    ib = mdiis_update(mdiis, ddcr, res, err);

    if ( verbose >= 2 )
      fprintf(stderr, "it %d, err %g -> %g, ib %d -> %d\n",
          it, (double) errp, (double) err, ibp, ib);
    if ( err < tol ) break;
    ibp = ib;
    errp = err;
  }
  /* use the best cr discovered so far */
  COPY1DARR(ddcr, mdiis->crbest, npt);
  /* update the corresponding ck, tr, tk */
  err = stepdd_picard(sphr, rho, ddcr, ddtr, ddck, ddtk,
      dtr, dck, dtk, tr, ck, tk, bphilr, fr, NULL, 0);
  if ( verbose || err > tol )
    fprintf(stderr, "iterdd_mdiis finished in %d steps, err %g\n", it, (double) err);

  mdiis_close(mdiis);
  return err;
}



#endif /* defined(IEMDIIS_H__) */

