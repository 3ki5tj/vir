#ifndef IEUTIL_H__
#define IEUTIL_H__



/* utility functions */



#ifndef PI
#define PI (xdouble) 3.1415926535897932384626433832795L
#endif

#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (int) (n); i_++) arr[i_] = 0; }

#define FREE1DARR(arr, n) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) arr[i_] = 0; \
  free(arr); }

#define COPY1DARR(dest, src, n) { int i_; \
  for ( i_ = 0; i_ < (n); i_++ ) dest[i_] = src[i_]; }

#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr, n1, n2) { \
  FREE1DARR(arr[0], (n1) * (n2)); \
  free(arr); }


__inline static xdouble integr(int n, xdouble *f, xdouble *w)
{
  xdouble y = 0;
  int i;

  for ( i = 0; i < n; i++ ) y += f[i] * w[i];
  return y;
}



__inline static xdouble integr2(int n, xdouble *f1, xdouble *f2, xdouble *w)
{
  xdouble y = 0;
  int i;

  for ( i = 0; i < n; i++ ) y += f1[i] * f2[i] * w[i];
  return y;
}



__inline static xdouble contactv(xdouble *f, int dm, xdouble B2)
{
  return B2 * (f[dm-1] + f[dm]) / 2;
}



/* compute t(k) from c(k) from the Ornstein-Zernike relation */
static void get_tk_oz(int l, int npt, xdouble **ck, xdouble **tk)
{
  int i, u;

  for ( i = 0; i < npt; i++ ) tk[l][i] = 0;
  for ( u = 0; u < l; u++ )
    for ( i = 0; i < npt; i++ )
      tk[l][i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]);
}



/* update the cavity function y(r)
 * for the hypernetted chain (HNC) approximation */
static void get_yr_hnc(int l, int nmax, int npt,
    xdouble **yr, xdouble *trl)
{
  int i, j, k, jl;
  xdouble *yro, powtr;

  xnew(yro, nmax - 1);
  /* y(r) = exp(t(r))
   * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
  /* yr = yr0 * exp(rho^l t_l)
   *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
   *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
   * y[l...n](r) are updated */
  for ( i = 0; i < npt; i++) {
    for (j = 0; j < nmax - 1; j++) yro[j] = yr[j][i];
    powtr = 1;
    for ( j = 1; j * l < nmax - 1; j++ ) {
      /* powtr = tl(r)^j/j! */
      powtr *= trl[i]/j;
      for ( jl = j * l, k = 0; k + jl < nmax - 1; k++ )
        /* yr_{k + jl} += yr0_k * tl^j/j! */
        yr[jl+k][i] += yro[k] * powtr;
    }
  }
  free(yro);
}



/* update the cavity function y(r)
 * for the Percus-Yevick approximation */
__inline static void get_yr_py2(int l, int npt, xdouble *yrl, xdouble **tr)
{
  int i, j;

  /* y(r) = 1 + t(r) + t(r)^2/2
   * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
  for ( i = 0; i < npt; i++) {
    yrl[i] = tr[l][i]; /* t(r) */
    for ( j = 1; j < l; j++ ) /* do t(r)^2 */
      yrl[i] += tr[j][i] * tr[l-j][i] / 2;
  }
}



/* correct the hypernetted chain approximation
 * `vc' is the trial correction function to y(r)
 * The corresponding correction to c(r) is (f(r) + 1) vc(r)
 * which is only effective for r > 1 */
__inline static xdouble get_corr1_hs(int l, int npt, int dm,
    xdouble *yr, xdouble *cr, xdouble *fr, xdouble *rDm1,
    xdouble B2, xdouble *vc, xdouble *Bc0, xdouble *Bv0, xdouble *eps)
{
  int i;
  xdouble dBc, dBv;

  /* Bc = - Int c(r) dr / (l + 2) */
  *Bc0 = dBc = 0;
  for ( i = 0; i < npt; i++ ) {
    *Bc0 += cr[i] * rDm1[i];
    dBc += vc[i] * (1 + fr[i]) * rDm1[i];
  }
  *Bc0 /= -(l + 2);
  dBc /= -(l + 2);
  *Bv0 = contactv(yr, dm, B2);
  dBv = contactv(vc, dm, B2);
  if ( l <= 1 ) {
    *eps = 0;
    return *Bc0;
  }
  *eps = -(*Bv0 - *Bc0) / (dBv - dBc);
  //printf("l %d, eps %g, Bv0 %g, Bc0 %g, dBv %g, dBc %g, B2 %g\n", l, (double) *eps, (double) *Bv0, (double) *Bc0, (double) dBv, (double) dBc, (double) B2);
  for ( i = 0; i < npt; i++ ) {
    vc[i] *= *eps;
    cr[i] += (fr[i] + 1) * vc[i];
  }
  return (*Bv0) + (*eps) * dBv;
}



/* correct the hypernetted chain approximation
 *  `vc' is the trial correction function to y(r) */
__inline static xdouble get_corr1(int l, int npt,
    xdouble *yr, xdouble *cr, xdouble *fr, xdouble *dfr,
    xdouble *rDm1, int dim, xdouble *vc,
    xdouble *Bc0, xdouble *Bv0, xdouble *eps)
{
  int i;
  xdouble dBc, dBv;

  *Bc0 = dBc = *Bv0 = dBv = 0;
  for ( i = 0; i < npt; i++ ) {
    *Bc0 += cr[i] * rDm1[i];
    dBc += vc[i] * (1 + fr[i]) * rDm1[i];
    *Bv0 += dfr[i] * yr[i] * rDm1[i];
    dBv += dfr[i] * vc[i] * rDm1[i];
  }
  *Bc0 /= -(l + 2);
  dBc /= -(l + 2);
  *Bv0 /= dim * 2;
  dBv /= dim * 2;
  if ( l <= 1 ) {
    *eps = 0;
    return *Bc0;
  }

  *eps = -(*Bv0 - *Bc0) / (dBv - dBc);
  for ( i = 0; i < npt; i++ ) {
    vc[i] *= *eps;
    cr[i] += vc[i] * (1 + fr[i]);
  }
  return *Bv0 + (*eps) * dBv;
}



/* Sum_{m = 3 to infinity} { rho^m [c(k)]^m }_{l+2} / m */
__inline static xdouble get_ksum(int l, int npt,
    xdouble **ck, xdouble *w, xdouble *s2)
{
  int i, j, k, m;
  xdouble *a, y, s = 0, z;

  if ( s2 != NULL ) *s2 = 0;
  MAKE1DARR(a, l + 1);
  for ( i = 0; i < npt; i++ ) {
    y = z = 0;
    for ( j = 0; j <= l; j++ )
      a[j] = ck[j][i];
    for ( m = 2; m <= l + 2; m++ ) { /* rho^m [c(k)]^m */
      /* update a[l + 2 - m], ..., a[0] */
      for ( j = l + 2 - m; j >= 0; j-- ) {
        /* a'[j] = Sum_{k = 0 to j} a[k] ck[j - k]
         * we start from large j to preserve small-index data */
        a[j] = a[j] * ck[0][i];
        for ( k = j - 1; k >= 0; k-- )
          a[j] += a[k] * ck[j - k][i];
      }
      if ( m == 2 ) continue;
      y += a[l + 2 - m] / m;
      z += a[l + 2 - m] * (m - 2) / (2*m);
    }
    s += w[i] * y;
    if ( s2 != NULL ) *s2 += w[i] * z;
  }
  FREE1DARR(a, l + 1);
  return s;
}



/* compute mu = Int (c - h t / 2) dr */
__inline static xdouble get_Bm_singer(int l, int npt,
    xdouble **cr, xdouble **tr, xdouble *w)
{
  int i, u;
  xdouble s = 0, x;

  for ( i = 0; i < npt; s += x * w[i], i++ )
    for ( x = cr[l][i], u = 0; u < l; u++ )
      x -= tr[l-u][i] * (cr[u][i] + tr[u][i]) / 2;
  return -s * (l+1)/(l+2);
}



/* compute -beta F = (1/2) Int (c - c t - t^2 / 2) dr (real part) */
__inline static xdouble get_Bh_singer(int l, int npt,
    xdouble **cr, xdouble **tr, xdouble *w)
{
  int i, u;
  xdouble s = 0, x;

  /* compute -2 beta F */
  for ( i = 0; i < npt; s += x * w[i], i++ )
    for ( x = cr[l][i], u = 0; u < l; u++ )
      x -= tr[l-u][i] * (cr[u][i] + tr[u][i] / 2);
  return -s * (l+1)/2;
}



/* d^2_rho beta P = Int (c + h t) dr */
__inline static xdouble get_Bx_py(int l, int npt,
    xdouble **cr, xdouble **tr, xdouble *w)
{
  int i, u;
  xdouble s = 0, x;

  for ( i = 0; i < npt; s += x * w[i], i++ )
    for ( x = cr[l][i], u = 0; u < l; u++ )
      x += tr[l-u][i] * (cr[u][i] + tr[u][i]);
  return -s/((l+1)*(l+2));
}



/* d^2_rho beta P = (-1/2) Int (c - t c) dr (real part) */
__inline static xdouble get_Bp_py(int l, int npt,
    xdouble **cr, xdouble **tr, xdouble *w)
{
  int i, u;
  xdouble s = 0, x;

  for ( i = 0; i < npt; i++ ) {
    for ( x = cr[l][i], u = 0; u < l; u++ )
      x -= tr[l-u][i] * cr[u][i];
    s += x * w[i];
  }
  return -s/2;
}



/* -Int h t dr / (l + 2) / l */
__inline static xdouble get_ht(int l, int npt,
    xdouble **cr, xdouble **tr, xdouble *w)
{
  int i, u;
  xdouble s = 0, x;

  for ( i = 0; i < npt; s += x * w[i], i++ )
    for ( x = 0, u = 0; u < l; u++ )
      x += tr[l-u][i] * (cr[u][i] + tr[u][i]);
  return -s/l/(l+2);
}



/* save the header for the virial file */
__inline static char *savevirhead(const char *fn, const char *title,
    int dim, int nmax, int doHNC, int mkcorr, int npt, xdouble rmax)
{
  FILE *fp;
  static char fndef[256];

  if ( fn == NULL ) {
    sprintf(fndef, "%sBn%s%sD%dn%dR%.0fM%d.dat",
        title ? title : "", doHNC ? "HNC" : "PY", mkcorr ? "c" : "",
        dim, nmax, (double) rmax, npt);
    fn = fndef;
  }
  xfopen(fp, fn, "w", return NULL);
  fprintf(fp, "# %s %s %d %.14f %d | n Bc Bv Bm [Bh Br | corr]\n",
      doHNC ? "HNC" : "PY", mkcorr ? "corr" : "",
      nmax, (double) rmax, npt);
  fclose(fp);
  return (char *) fn;
}



/* print a virial coefficient */
__inline static int printB(const char *name, int dim, int n,
    xdouble B, xdouble B2p, xdouble volp, const char *ending)
{
  double x;
  int px = 0;

  if ( B == 0 ) return px;
  printf("%s(%3d) = %15.8e", name, n, (double) B);
  if ( B2p != 0 ) {
    printf(" (%15.8e", (double) (B/B2p));
    x = (double) (B/volp);
    px = ( fabs(x) < 10000 && dim % 2 == 1);
    if ( !px ) printf(")");
    else printf(",%+12.6f)", x);
  }
  printf("%s", ending);
  return px;
}



/* save virial coefficients */
__inline static int savevir(const char *fn, int dim, int n,
    xdouble Bc, xdouble Bv, xdouble Bm, xdouble Bh, xdouble Br,
    xdouble B2p, int mkcorr, xdouble fcorr)
{
  FILE *fp;
  xdouble volp;

  /* print the result on screen */
  volp = B2p / pow_si(2, (n-1)*(dim-1));
  printB("Bc", dim, n, Bc, B2p, volp, ", ");
  printB("Bv", dim, n, Bv, B2p, volp, ", ");
  /* when making corrections, Bm is the corrected value */
  printB("Bm", dim, n, Bm, B2p, volp, "");
  if ( mkcorr ) {
    printf(", %9.6f\n", (double) fcorr);
  } else { /* the following are useless when making corrections */
    printf("\n");
    if ( Bh != 0 || Br != 0 ) {
      printB("Bh", dim, n, Bh, B2p, volp, ", ");
      printB("Br", dim, n, Br, B2p, volp, "\n");
    }
  }

  if (fn != NULL) {
    xfopen(fp, fn, "a", return -1);
    /* normalize the virial coefficients */
    if ( B2p != 0 ) {
      Bc /= B2p;
      Bv /= B2p;
      Bm /= B2p;
      Bh /= B2p;
      Br /= B2p;
    }
    if ( mkcorr ) {
      fprintf(fp, "%4d%+24.14e%+24.14e%+24.14e %+18.14f\n",
          n, (double) Bc, (double) Bv, (double) Bm, (double) fcorr);
    } else {
      fprintf(fp, "%4d%+24.14e%24.14e%24.14e%24.14e%24.14e\n",
          n, (double) Bc, (double) Bv,
          (double) Bm, (double) Bh, (double) Br);
    }
    fclose(fp);
  }

  return 0;
}



/* save c(r) or t(r) file */
__inline static int savecrtr(const char *fn, int l, int npt,
    xdouble *ri, xdouble *cr, xdouble *tr, xdouble *vc, xdouble **yr)
{
  FILE *fp;
  int i, j;

  if (fn == NULL) return -1;
  xfopen(fp, fn, (l == 1) ? "w" : "a", return -1);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%10.7f %20.12f %20.12f %d",
        (double) ri[i], (double) cr[i], (double) tr[i], l);
    if ( vc != NULL )
      fprintf(fp, " %20.12f", (double) vc[i]);
    if ( yr != NULL )
      for ( j = 1; j <= l; j++ )
        fprintf(fp, " %20.12f", (double) yr[j][i]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  fclose(fp);
  return 0;
}



#endif
