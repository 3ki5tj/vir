#ifndef IEUTIL_H__
#define IEUTIL_H__



/* utility functions */



/* in case xdouble hasn't been included, include it here */
#include "xdouble.h"



#if HAVEF128
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wformat-extra-args"
#endif



#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (int) (n); i_++) arr[i_] = 0; }

#define FREE1DARR(arr, n) if ( (arr) != NULL ) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) (arr)[i_] = 0; \
  free(arr); }

#define COPY1DARR(dest, src, n) { int i_; \
  for ( i_ = 0; i_ < (n); i_++ ) dest[i_] = src[i_]; }

#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr, n1, n2) if ( (arr) != NULL ) { \
  FREE1DARR((arr)[0], (n1) * (n2)); \
  free(arr); }



__inline static xdouble adjustrmax(xdouble rmax, int npt,
    xdouble *dr, int *dm, int ffttype)
{
  /* fix dr such that r = 1 at the middle of bins dm - 1 and dm */
  if ( ffttype ) { /* r(i) = dr * (i + .5) */
    *dm = (int) (npt / rmax + 1e-8); /* number of bins in the core */
    *dr = (xdouble) 1 / *dm;
  } else {
    *dm = (int) (npt / rmax + .5 + 1e-8);
    *dr = (xdouble) 2 / (*dm*2 - 1);
  }
  rmax = *dr * npt;
  printf("%d bins, %d within the hard core (dr %g), rmax %g\n",
      npt, *dm, (double) *dr, (double) rmax);
  return rmax;
}



#ifdef SIMPSON
/* the Simpson formula does not improve the precision too much */
#define INT_GETW(i) ((xdouble) (i % 2 ? 4 : i == 0 ? 1 : 2) / 3)
#else
#define INT_GETW(i) 1
#endif


__inline static xdouble integr(int n, xdouble *f, xdouble *w)
{
  xdouble y = 0;
  int i;

  for ( i = 0; i < n; i++ )
    y += f[i] * w[i] * INT_GETW(i);
  return y;
}



__inline static xdouble integr2(int n, xdouble *f1, xdouble *f2, xdouble *w)
{
  xdouble y = 0;
  int i;

  for ( i = 0; i < n; i++ )
    y += f1[i] * f2[i] * w[i] * INT_GETW(i);
  return y;
}



__inline static xdouble integre(int n, xdouble *f, xdouble *fr, xdouble *w)
{
  xdouble y = 0;
  int i;

  for ( i = 0; i < n; i++ )
    y += f[i] * (fr[i] + 1) * w[i] * INT_GETW(i);
  return y;
}



__inline static xdouble contactv(xdouble *f, int dm, xdouble B2)
{
  return B2 * (f[dm-1] + f[dm]) / 2;
}



__inline static xdouble get_Bv_hnc_hs(xdouble *yrl, xdouble *trl,
    xdouble a0, xdouble q, int dm, xdouble B2)
{
  xdouble y = contactv(yrl, dm, B2);
  if (yrl != trl) /* to avoid PY case */
    y += (1 - a0*q) * contactv(trl, dm, B2);
  return y;
}



__inline static xdouble get_Bv_hnc(int npt, xdouble *yrl, xdouble *trl,
    xdouble a0, xdouble q, xdouble *rdfr, xdouble *rDm1, int dim)
{
  xdouble y = integr2(npt, yrl, rdfr, rDm1);
  if (yrl != trl)
    y += (1 - a0*q) * integr2(npt, trl, rdfr, rDm1);
  return y / (dim * 2);
}



/* compute t(k) from c(k) from the Ornstein-Zernike relation */
__inline static void get_tk_oz(int l, int npt, xdouble **ck, xdouble **tk)
{
  int i, u;

  for ( i = 0; i < npt; i++ ) tk[l][i] = 0;
  for ( u = 0; u < l; u++ )
    for ( i = 0; i < npt; i++ )
      tk[l][i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]);
}



#define get_yr_hnc(l, nmax, npt, yr, trl) \
  get_yr_hncx(l, nmax, npt, yr, trl, 1, NULL)

/* update the cavity function y(r)
 * for the generalized hypernetted chain (HNC) approximation */
__inline static void get_yr_hncx(int l, int nmax, int npt,
    xdouble **yr, xdouble *trl, xdouble q, xdouble *swr)
{
  int i, j, k, jl;
  xdouble *yro, powtr;

  xnew(yro, nmax - 1);
  /* y(r) = a0 exp[q sw(r) t(r) ] / sw(r)
   * where sw(r) = 1 - exp(-alpha r)
   * c(r) = (1 + f(r)) (y(r) + (1-a0/sw(r)) + (1-a0*q)*t(r)) - t(r) - 1 */
  /* yr = yr0 * exp(q sw rho^l t_l)
   *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
   *    * (1 + q sw t_l rho^l + q^2 sw^2 t_{2l} rho^(2l)/2! + ... )
   * y[l...n](r) are updated */
  for ( i = 0; i < npt; i++) {
    for (j = 0; j < nmax - 1; j++) yro[j] = yr[j][i];
    powtr = 1;
    for ( j = 1; j * l < nmax - 1; j++ ) {
      /* powtr = q^j tl(r)^j/j! */
      powtr *= q * trl[i] / j;
      if ( swr != NULL ) powtr *= swr[i];
      for ( jl = j * l, k = 0; k + jl < nmax - 1; k++ )
        /* yr_{k + jl} += yr0_k * q^j tl^j/j! */
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
 * which is only effective for r > 1
 * i.e., vc(r < 1) is unused */
__inline static xdouble get_corr1_hs(int l, int npt, int dm,
    xdouble *cr, xdouble *fr, xdouble *rDm1,
    xdouble B2, xdouble *vc, xdouble *Bc0, xdouble *Bv0, xdouble *eps)
{
  int i;
  xdouble dBc, dBv;

  //*Bc0 = -integr(npt, cr, rDm1) / (l + 2);
  dBc = -integre(npt, vc, fr, rDm1) / (l + 2);
  //*Bv0 = contactv(yr, dm, B2);
  dBv = contactv(vc, dm, B2);
  if ( l <= 1 ) {
    *eps = 0;
    return *Bv0;
  }
  *eps = -(*Bv0 - *Bc0) / (dBv - dBc);
  for ( i = 0; i < npt; i++ ) {
    vc[i] *= *eps;
    cr[i] += (fr[i] + 1) * vc[i];
  }
  return (*Bv0) + (*eps) * dBv;
}



/* switch between the continuous and discontinuous case */
#define get_corr1x(l, npt, dm, cr, fr, rdfr, rDm1, dim, B2, vc, Bc0, Bv0, fcorr) \
  (rdfr == NULL \
    ? get_corr1_hs(l, npt, dm, cr, fr, rDm1, B2, vc, Bc0, Bv0, fcorr) \
    : get_corr1(l, npt, cr, fr, rdfr, rDm1, dim, vc, Bc0, Bv0, fcorr))



/* correct the hypernetted chain approximation
 *  `vc' is the trial correction function to y(r) */
__inline static xdouble get_corr1(int l, int npt,
    xdouble *cr, xdouble *fr, xdouble *rdfr,
    xdouble *rDm1, int dim, xdouble *vc,
    xdouble *Bc0, xdouble *Bv0, xdouble *eps)
{
  int i;
  xdouble dBc, dBv;

  //*Bc0 = -integr(npt, cr, rDm1) / (l + 2);
  dBc = -integre(npt, vc, fr, rDm1) / (l + 2);
  //*Bv0 = integr2(npt, yr, rdfr, rDm1) / (dim * 2);
  dBv = integr2(npt, vc, rdfr, rDm1) / (dim * 2);
  if ( l <= 1 ) {
    *eps = 0;
    return *Bc0;
  }

  *eps = -(*Bv0 - *Bc0) / (dBv - dBc);
  *eps = 1;
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



#define savevirhead(fn, title, dim, l0, nmax, dohnc, mkcorr, npt, rmax, inittime) \
  savevirheadx(fn, title, dim, l0, nmax, dohnc, mkcorr, npt, rmax, inittime, \
      1, 1, -1, 0)

/* save the header for the virial file */
__inline static char *savevirheadx(const char *fn, const char *title,
    int dim, int l0, int nmax, int dohnc, int mkcorr, int npt,
    xdouble rmax, clock_t inittime,
    xdouble hncamp, xdouble hncq, xdouble hncalpha, xdouble shift)
{
  FILE *fp;
  static char fndef[512];

  if ( fn == NULL ) {
    char shncamp[80] = "", shncq[80] = "", shncalpha[80] = "";
    char sshift[80] = "";
    if ( FABS(hncamp - 1) > 1e-6 )
      sprintf(shncamp, "a%g", (double) hncamp);
    if ( FABS(hncq - 1) > 1e-6 )
      sprintf(shncq, "q%g", (double) hncq);
    if ( hncalpha >= 0 )
      sprintf(shncalpha, "A%g", (double) hncalpha);
    if ( FABS(shift) > 1e-6 )
      sprintf(sshift, "c%g", (double) shift);
    sprintf(fndef, "%sBn%s%sD%dn%dR%.0fM%d%s%s%s%s%s.dat",
        title ? title : "", dohnc ? "HNC" : "PY", mkcorr ? "c" : "",
        dim, nmax, (double) rmax, npt,
        shncamp, shncq, shncalpha, sshift,
        STRPREC);
    fn = fndef;
  }
  xfopen(fp, fn, (l0 == 1) ? "w" : "a", return NULL);
  fprintf(fp, "# %s %s %d %.14f %d %s | n Bc Bv Bm [Bh Br | corr] | %.3fs\n",
      dohnc ? "HNC" : "PY", mkcorr ? "corr" : "",
      nmax, (double) rmax, npt, STRPREC,
      (double) inittime / CLOCKS_PER_SEC);
  fclose(fp);
  return (char *) fn;
}



/* print a virial coefficient */
__inline static int printB(const char *name, int dim, int n,
    xdouble B, xdouble B2p, xdouble B2q,
    xdouble volp, xdouble volq, const char *ending)
{
  xdouble x;
  int px = 0;

  if ( B == 0 ) return px;
  printf("%s(%3d) = %15.7" XDBLPRNF "e", name, n, B);
  if ( B2p != 0 ) {
    printf(" (%15.8" XDBLPRNF "e", B/B2p/B2q);
    x = B/volp/volq;
    px = ( FABS(x) < 10000 && dim % 2 == 1);
    if ( !px ) printf(")");
    else printf(",%+12.6" XDBLPRNF "f)", x);
  }
  printf("%s", ending);
  return px;
}



/* save a virial coefficient to file */
__inline static void saveB(FILE *fp, xdouble B, xdouble B2p, xdouble B2q)
{
  if ( B == 0 ) fprintf(fp, " 0");
  else fprintf(fp, " " XDBLPRNE, (B2p != 0 ? B/B2p/B2q : B));
}



/* save virial coefficients */
__inline static int savevir(const char *fn, int dim, int n,
    xdouble Bc, xdouble Bv, xdouble Bm, xdouble Bh, xdouble Br,
    xdouble B2, int mkcorr, xdouble fcorr)
{
  FILE *fp;
  int np = n - 1, nq = 0;
  xdouble vol, volp, volq = 1, B2p, B2q = 1;

  /* print the result on screen */
  B2p = pow_si(B2, np);
  if ( B2p == 0 ) { /* underflow */
    nq = np / 2;
    np -= nq;
    B2p = pow_si(B2, np);
    B2q = pow_si(B2, nq);
  }
  vol = B2/pow_si(2, dim-1);
  volp = pow_si(vol, np);
  volq = pow_si(vol, nq);
  printB("Bc", dim, n, Bc, B2p, B2q, volp, volq, ", ");
  printB("Bv", dim, n, Bv, B2p, B2q, volp, volq, ", ");
  /* when making corrections, Bm is the corrected value */
  printB("Bm", dim, n, Bm, B2p, B2q, volp, volq, "");
  if ( mkcorr ) {
    printf(", %9.6f\n", (double) fcorr);
  } else { /* the following are useless when making corrections */
    printf("\n");
    if ( Bh != 0 || Br != 0 ) {
      printB("Bh", dim, n, Bh, B2p, B2q, volp, volq, ", ");
      printB("Br", dim, n, Br, B2p, B2q, volp, volq, "\n");
    }
  }

  if (fn != NULL) {
    xfopen(fp, fn, "a", return -1);
    fprintf(fp, "%4d", n);
    saveB(fp, Bc, B2p, B2q);
    saveB(fp, Bv, B2p, B2q);
    saveB(fp, Bm, B2p, B2q);
    if ( mkcorr ) {
      fprintf(fp, " " XDBLPRNE "\n", fcorr);
    } else {
      saveB(fp, Bh, B2p, B2q);
      saveB(fp, Br, B2p, B2q);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  return 0;
}



/* save head */
__inline static int savevirtail(const char *fn, clock_t endtime)
{
  FILE *fp;

  if (fn == NULL) return -1;
  xfopen(fp, fn, "a", return -1);
  fprintf(fp, "# %.3fs\n", (double) endtime / CLOCKS_PER_SEC);
  fclose(fp);
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
    fprintf(fp, XDBLPRNE " " XDBLPRNE " " XDBLPRNE " %d",
        ri[i], cr[i], tr[i], l);
    if ( vc != NULL )
      fprintf(fp, " " XDBLPRNE, vc[i]);
    if ( yr != NULL )
      for ( j = 1; j <= l; j++ )
        fprintf(fp, " " XDBLPRNE, yr[j][i]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  fclose(fp);
  return 0;
}



/* snapshot variables */
char fnck[80], fntk[80], fncr[80], fntr[80], fnyrstem[80], fnyr[80];

#define SNAPSHOT_READARRTXT(fp, fn, arr, n) { int i_; \
  for ( i_ = 0; i_ < n; i_++ ) \
    fscanf(fp, "%" XDBLSCNF "f", (xdouble *)(arr + i_)); }

#define SNAPSHOT_WRITEARRTXT(fp, fn, arr, n) { int i_; \
  for ( i_ = 0; i_ < n; i_++ ) \
    fprintf(fp, "%22.14" XDBLPRNF "e ", *((xdouble *)(arr)+i_)); \
  fprintf(fp, "\n"); }

#define SNAPSHOT_READARRBIN(fp, fn, arr, n) \
  die_if ( fread(arr, sizeof(xdouble), n, fp) != (size_t) n, \
    "cannot read %ld elements from %s\n", (long) n, fn)

#define SNAPSHOT_WRITEARRBIN(fp, fn, arr, n) \
  die_if ( fwrite(arr, sizeof(xdouble), n, fp) != (size_t) n, \
    "cannot write %ld elements to %s\n", (long) n, fn)

#ifdef SNAPSHOT_TEXTIO
int snapshot_textio = 1;
#define SNAPSHOT_READARR    SNAPSHOT_READARRTXT
#define SNAPSHOT_WRITEARR   SNAPSHOT_WRITEARRTXT
#else
int snapshot_textio = 0;
#define SNAPSHOT_READARR    SNAPSHOT_READARRBIN
#define SNAPSHOT_WRITEARR   SNAPSHOT_WRITEARRBIN
#endif



/* try to load array from file `fn' */
__inline static int snapshot_loadf(int *l, int *l0, int nmax,
    const char *fn, int npt, xdouble **arr, int offset, xdouble *arrl)
{
  FILE *fp;

  if ( (fp = fopen(fn, "r+b")) == NULL ) {
    fprintf(stderr, "cannot load %s\n", fn);
    return -1;
  }

  if ( snapshot_textio ) { /* text file */
    for ( *l = 0; ; ) {
      if (fgetc(fp) == '\n') (*l)++;
      if (feof(fp)) break;
    }
  } else { /* binary file */
    long fpos, arrsz = npt * (long) sizeof(xdouble);
    fseek(fp, 0, SEEK_END);
    fpos = ftell(fp);
    die_if ( fpos % arrsz != 0,
      "%s is corrupted %ld / %ld\n", fn, fpos, arrsz);
    /* check the size */
    *l = (int) (fpos / arrsz);
    if ( *l0 < 0 ) {
      *l0 = *l;
    } else { /* check against the known value */
      die_if ( *l0 != *l,
        "%s has %d arrays, instead of %d\n", fn, *l, *l0);
    }
  }

  if ( arr != NULL && *l > 0 ) { /* read the 2D array */
    int l1 = (*l >= nmax - 1 - offset) ? nmax - 1 - offset : *l;
    rewind(fp);
    SNAPSHOT_READARR(fp, fn, arr[offset], npt*l1);
  }

  if ( arrl != NULL && *l > 0 ) { /* read the last 1D array */
    if ( snapshot_textio ) { /* read text file */
      int j;
      rewind(fp);
      for ( j = 0; j < *l; j++ )
        SNAPSHOT_READARR(fp, fn, arrl, npt);
    } else { /* read binary file */
      fseek(fp, (*l-1)*npt*sizeof(xdouble), SEEK_SET);
      SNAPSHOT_READARR(fp, fn, arrl, npt);
    }
    /*printf("%s, %s l %d\n", fn, #arrl, *l);*/
  }
  fclose(fp);
  return 0;
}



__inline static int snapshot_open(int dim, int nmax, xdouble rmax,
    int dohnc, int mkcorr, int ring, int singer, int npt,
    xdouble **ck, xdouble **tk, xdouble **cr, xdouble **tr,
    xdouble *crl, xdouble *trl, xdouble **yr)
{

#define SNAPSHOT_MKSTEM(fn, arr) \
  sprintf(fn, "snapshotD%d%s%s%s%sR%.3fM%d%s_%s", \
      dim, dohnc ? "HNC" : "PY", mkcorr ? "c" : "", \
      ring ? "r" : "", singer ? "s" : "", \
      (double) rmax, npt, STRPREC, #arr)

#define SNAPSHOT_OPENF(l, l0, fn, arr, offset, arrl) { \
  SNAPSHOT_MKSTEM(fn, arr); strcat(fn, ".dat"); \
  snapshot_loadf(&l, &l0, nmax, fn, npt, arr, offset, arrl); }

  int l = 0, l0 = -1, n1 = 0, n2 = nmax - 1;

  SNAPSHOT_OPENF(l, l0, fnck, ck, 0, NULL);
  SNAPSHOT_OPENF(l, l0, fntk, tk, 1, NULL);
  SNAPSHOT_OPENF(l, l0, fncr, cr, 1, crl);
  SNAPSHOT_OPENF(l, l0, fntr, tr, 1, trl);
  if ( yr != NULL ) {
    SNAPSHOT_MKSTEM(fnyrstem, yr);
    if ( l > 0 ) {
      sprintf(fnyr, "%sl%d.dat", fnyrstem, l);
      if ( snapshot_loadf(&n1, &n2, nmax, fnyr, npt, yr, 0, NULL) != 0 ) {
        fprintf(stderr, "cannot load yr(%d) of %d arrays from %s\n",
            l, n2, fnyr);

#if 1
        /* try the old file name, NOTE: only for compatibility */
        sprintf(fnyr, "%s.dat", fnyrstem);
        n1 = 0; n2 = nmax - 1;
        if ( snapshot_loadf(&n1, &n2, nmax, fnyr, npt, yr, 0, NULL) != 0 ) {
          fprintf(stderr, "cannot load yr (old style) of %d arrays from %s\n", n2, fnyr);
          exit(-1);
        }
#else
        /* currently we simply die if fnyr doesn't exists
         * however, yr can be reconstructed, as in a rerun */
        exit(-1);
#endif

      }
    } /* l > 0 */
  }
  printf("snapshot initial l %d/%d\n", l+1, n2);
  return l + 1;
}



/* append an array `arr' of size `npt' to file `fn' */
__inline static int snapshot_appendf(const char *fn, xdouble *arr, int npt)
{
  if ( arr != NULL ) {
    FILE *fp;

    xfopen(fp, fn, "ab", exit(-1));
    SNAPSHOT_WRITEARR(fp, fn, arr, npt);
    fclose(fp);
  }
  return 0;
}



__inline static void snapshot_take(int l, int npt,
    xdouble *ckl, xdouble *tkl, xdouble *crl, xdouble *trl,
    int nmax, xdouble **yr)
{
  snapshot_appendf(fnck, ckl, npt);
  snapshot_appendf(fntk, tkl, npt);
  snapshot_appendf(fncr, crl, npt);
  snapshot_appendf(fntr, trl, npt);
  if ( yr != NULL ) {
    FILE *fp;
    int i;

    sprintf(fnyr, "%sl%d.dat", fnyrstem, l);
    xfopen(fp, fnyr, "wb", exit(-1));
    for ( i = 0; i < nmax - 1; i++)
      SNAPSHOT_WRITEARR(fp, fnyr, yr[i], npt);
    fclose(fp);
  }
}



#if HAVEF128
#pragma GCC diagnostic pop
#endif



#endif
