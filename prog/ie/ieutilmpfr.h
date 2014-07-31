#ifndef IEUTIL_MPFR_H
#define IEUTIL_MPFR_H


#define INIT_(a) mpfr_init(a)

#define CLEAR_(a) mpfr_clear(a)

/* y = x */
#define SET_(y, x) mpfr_set(y, x, MPFR_RNDN)

/* x = i */
#define SET_SI_(x, i) mpfr_set_si(x, i, MPFR_RNDN)

/* x = s */
#define SET_STR_(x, s) mpfr_set_str(x, s, 10, MPFR_RNDN)

/* z = x + y */
#define ADD_(z, x, y) mpfr_add(z, x, y, MPFR_RNDN)

/* y += x */
#define ADD_X_(y, x) ADD_(y, y, x)

/* y = x + i */
#define ADD_SI_(y, x, i) mpfr_add_si(y, x, i, MPFR_RNDN)

/* y += i */
#define ADD_SI_X_(y, i) ADD_SI_(y, y, i)

/* z = x - y */
#define SUB_(z, x, y) mpfr_sub(z, x, y, MPFR_RNDN)

/* y -= x */
#define SUB_X_(y, x) SUB_(y, y, x)

/* y = x - i */
#define SUB_SI_(y, x, i) mpfr_sub_si(y, x, i, MPFR_RNDN)

/* y -= i */
#define SUB_SI_X_(y, i) SUB_SI_(y, y, i)

/* y = -x */
#define NEG_(y, x) mpfr_neg(y, x, MPFR_RNDN)

/* z = x * y */
#define MUL_(z, x, y) mpfr_mul(z, x, y, MPFR_RNDN)

/* y *= x */
#define MUL_X_(y, x) MUL_(y, y, x)

/* y = x * i */
#define MUL_SI_(y, x, i) mpfr_mul_si(y, x, i, MPFR_RNDN)

/* y *= i */
#define MUL_SI_X_(y, i) MUL_SI_(y, y, i)

/* y = a / b */
#define DIV_(y, a, b) mpfr_div(y, a, b, MPFR_RNDN)

/* y /= x */
#define DIV_X_(y, x) DIV_(y, y, x)

/* y = x / i */
#define DIV_SI_(y, x, i) mpfr_div_si(y, x, i, MPFR_RNDN)

/* y /= i */
#define DIV_SI_X_(y, i) DIV_SI_(y, y, i)

/* y = i / x */
#define SI_DIV_(y, i, x) mpfr_si_div(y, i, x, MPFR_RNDN)

/* y = i / y */
#define SI_DIV_X_(y, i) SI_DIV_(y, i, y)

/* y = a * b + x */
#define FMA_(y, a, b, x) mpfr_fma(y, a, b, x, MPFR_RNDN)

/* y += a * b */
#define FMA_X_(y, a, b) FMA_(y, a, b, y)

/* y = a * b - x */
#define FMS_(y, a, b, x) mpfr_fms(y, a, b, x, MPFR_RNDN)

/* y = x * x */
#define SQR_(y, x) mpfr_sqr(y, x, MPFR_RNDN)

/* x = x * x */
#define SQR_X_(x) SQR_(x, x)

/* xabs = |x| */
#define ABS_(xabs, x) mpfr_abs(xabs, x, MPFR_RNDN)

/* y = x ^ i */
#define POW_SI_(y, x, i) mpfr_pow_si(y, x, i, MPFR_RNDN)

/* x = x ^ i */
#define POW_SI_X_(x, i) POW_SI_(x, x, i)

#define CONST_PI_(pi) mpfr_const_pi(pi, MPFR_RNDN)

#define GET_D_(x) mpfr_get_d(x, MPFR_RNDN)

/* y = a * b * c */
#define MUL3_(y, a, b, c) { \
  mpfr_mul(y, a, b, MPFR_RNDN); mpfr_mul(y, y, c, MPFR_RNDN); }

/* y *= a * b */
#define MUL3_X_(y, a, b) MUL3_(y, y, a, b)

/* y = (a + b) * c */
#define FAM_(y, a, b, c) { \
  mpfr_add(y, a, b, MPFR_RNDN); mpfr_mul(y, y, c, MPFR_RNDN); }

/* sign of a - b */
#define CMP_(a, b) mpfr_cmp(a, b)

/* sign of x - i */
#define CMP_SI_(x, i) mpfr_cmp_si(x, i)



#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (int) (n); i_++) { \
    INIT_(arr[i_]); \
    SET_SI_(arr[i_], 0); } }

#define FREE1DARR(arr, n) if ((arr) != NULL) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) \
    CLEAR_((arr)[i_]); \
    free(arr); }

#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr, n1, n2) if ((arr) != NULL) { \
  FREE1DARR((arr)[0], (n1) * (n2)); \
  free(arr); }

#define COPY1DARR(a, b, n) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) \
    SET_(a[i_], b[i_]); }


__inline static double adjustrmax(const char *srmax, int npt,
    mpfr_t dr, int *dm, int ffttype)
{
  double rmax = atof(srmax);

  /* fix dr such that r = 1 at the middle of bins dm - 1 and dm */
  if ( ffttype ) { /* r(i) = dr * (i + .5) */
    *dm = (int) (npt / rmax + 1e-8); /* number of bins in the core */
    SET_SI_(dr, 1);
    DIV_SI_X_(dr, *dm);
  } else {
    *dm = (int) (npt / rmax + .5 + 1e-8);
    SET_SI_(dr, 2);
    DIV_SI_X_(dr, *dm * 2 - 1);
  }
  rmax = GET_D_(dr) * npt;
  printf("%d bins, %d within the hard core (dr %g), rmax %g\n",
      npt, *dm, GET_D_(dr), rmax);
  return rmax;
}



__inline static void integr(mpfr_t y, int n, mpfr_t *f, mpfr_t *w)
{
  int i;

  SET_SI_(y, 0);
  for ( i = 0; i < n; i++ )
    FMA_X_(y, f[i], w[i]); /* y += f[i] * w[i] */
}



__inline static void integr2(mpfr_t y, int n,
    mpfr_t *f, mpfr_t *f2, mpfr_t *w)
{
  int i;
  mpfr_t x;

  INIT_(x);
  SET_SI_(y, 0);
  for ( i = 0; i < n; i++ ) {
    MUL_(x, f[i], f2[i]);
    FMA_X_(y, x, w[i]); /* y += x * w[i] */
  }
  CLEAR_(x);
}



/* B2 * (f[dm-1] + f[dm]) / 2 */
__inline static void contactv(mpfr_t y, mpfr_t *f, int dm, mpfr_t B2)
{
  FAM_(y, f[dm-1], f[dm], B2);
  DIV_SI_X_(y, 2);
}



/* get the value at zero separation */
__inline static void get_zerosep(mpfr_t z, mpfr_t *y, mpfr_t *x)
{
  mpfr_t a, b, c;
  INIT_(a);
  INIT_(b);
  INIT_(c);
  // return (y[0]*x[1] - y[1]*x[0])/(x[1] - x[0]);
  MUL_(a, y[0], x[1]);
  MUL_(b, y[1], x[0]);
  SUB_(z, a, b);
  SUB_(c, x[1], x[0]);
  DIV_X_(z, c);
  CLEAR_(a);
  CLEAR_(b);
  CLEAR_(c);
}



/* return the virial coefficient of order l + 1 */
__inline static void update_lnyr0(mpfr_t B, int l, mpfr_t *yr0, mpfr_t *lnyr0)
{
  int j;
  mpfr_t a;

  INIT_(a);
  MUL_SI_(lnyr0[l], yr0[l], l); /* lnyr0[l] = l * yr0[l]; */
  for ( j = 1; j < l; j++ ) {
    /* lnyr0[l] -= j * lnyr0[j] * yr0[l - j]; */
    MUL_(a, lnyr0[j], yr0[l-j]);
    MUL_SI_X_(a, j);
    SUB_X_(lnyr0[l], a);
  }
  DIV_SI_(B, lnyr0[l], l + 1);
  DIV_SI_X_(lnyr0[l], l);
  CLEAR_(a);
}



/* compute t(k) from c(k) from the Ornstein-Zernike relation */
__inline static void get_tk_oz(int l, int npt, mpfr_t **ck, mpfr_t *tkl)
{
  int i, u, v;
  mpfr_t x, *a;

  INIT_(x);
  MAKE1DARR(a, l + 1);
  for ( i = 0; i < npt; i++ ) {
    for (v = 1; v <= l; v++) {
      SET_SI_(a[v], 0);
      for ( u = 0; u < v; u++ ) {
        /* a[v] += ck[v-1-u][i] * (ck[u][i] + a[u]) */
        ADD_(x, ck[u][i], a[u]);
        FMA_X_(a[v], ck[v-1-u][i], x);
      }
    }
    SET_(tkl[i], a[l]);
  }
  FREE1DARR(a, l + 1);
  CLEAR_(x);
}



/* compute t(k) from c(k) from the Ornstein-Zernike relation */
__inline static void get_tk_oz_fast(int l, int npt, mpfr_t **ck, mpfr_t **tk)
{
  int i, u;
  mpfr_t x;

  INIT_(x);
  for ( i = 0; i < npt; i++ )
    SET_SI_(tk[l][i], 0); /* tl[i] = 0; */
  for ( u = 0; u < l; u++ )
    for ( i = 0; i < npt; i++ ) {
      /* tl[i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]); */
      ADD_(x, ck[u][i], tk[u][i]);
      FMA_X_(tk[l][i], ck[l-1-u][i], x);
    }
  CLEAR_(x);
}



/* compute the cavity function y(r)
 * from the hypernetted chain (HNC) approximation */
static void get_yr_hnc(int l, int npt, mpfr_t *yrl, mpfr_t **tr)
{
  int i, u, v;
  mpfr_t x, *y;

  INIT_(x);
  MAKE1DARR(y, l + 1);
  /* y(r) = exp(t(r))
   * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
  for ( i = 0; i < npt; i++) {
    SET_SI_(y[0], 1);
    for ( v = 1; v <= l; v++ ) { /* compute y[v] */
      SET_SI_(y[v], 0);
      for ( u = 1; u <= v; u++ ) {
        /* y[v] += u * tr[u][i] * y[v - u]; */
        MUL_(x, tr[u][i], y[v - u]);
        MUL_SI_X_(x, u);
        ADD_X_(y[v], x);
      }
      /* y[v] /= v; */
      DIV_SI_X_(y[v], v);
    }
    SET_(yrl[i], y[l]);
  }
  CLEAR_(x);
  FREE1DARR(y, l + 1);
}



/* compute the cavity function y(r)
 * from the hypernetted chain (HNC) approximation */
static void get_yr_hnc_fast(int l, int npt, mpfr_t *yrl,
    mpfr_t **yr, mpfr_t **tr)
{
  int i, u;
  mpfr_t x;

  INIT_(x);
  /* y(r) = exp(t(r))
   * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
  for ( i = 0; i < npt; i++) {
    SET_SI_(yr[l][i], 0);
    for ( u = 1; u <= l; u++ ) {
      /* yr[l][i] += u * tr[u][i] * yr[l - u][i]; */
      MUL_(x, tr[u][i], yr[l - u][i]);
      MUL_SI_X_(x, u);
      ADD_X_(yr[l][i], x);
    }
    /* y[l] /= l; */
    DIV_SI_X_(yr[l][i], l);
    SET_(yrl[i], yr[l][i]);
  }
  CLEAR_(x);
}



/* correct the hypernetted chain approximation
 *  `vc' is the trial correction function to y(r) */
__inline static void get_corr1_hs(mpfr_t B, int l, int npt, int dm,
    mpfr_t *ylr, mpfr_t *crl, mpfr_t *fr, mpfr_t *rDm1,
    mpfr_t B2, mpfr_t *vc, mpfr_t Bc0, mpfr_t Bv0, mpfr_t eps)
{
  int i;
  mpfr_t dBc, dBv, x;

  INIT_(dBc);
  INIT_(dBv);
  INIT_(x);

  SET_SI_(Bc0, 0);
  SET_SI_(dBc, 0);
  for ( i = 0; i < npt; i++ ) {
    /* Bc0 += crl[i] * rDm1[i]; */
    FMA_X_(Bc0, crl[i], rDm1[i]);

    /* dBc += vc[i] * (1 + fr[i]) * rDm1[i]; */
    ADD_SI_(x, fr[i], 1);
    MUL3_X_(x, vc[i], rDm1[i]);
    ADD_X_(dBc, x);
  }
  DIV_SI_X_(Bc0, -(l+2)); /* Bc0 /= -(l+2); */
  DIV_SI_X_(dBc, -(l+2)); /* dBc /= -(l+2); */

  /* Bv0 = B2 * (ylr[dm] + ylr[dm-1])/2; */
  contactv(Bv0, ylr, dm, B2);
  /* dBv = B2 * (vc[dm] + vc[dm-1])/2; */
  contactv(dBv, vc, dm, B2);
  if ( l <= 1 ) {
    SET_(B, Bv0);
    SET_SI_(eps, 0);
    return;
  }

  /* eps = -(Bv0 - Bc0) / (dBv - dBc); */
  SUB_(eps, Bc0, Bv0);
  SUB_(x, dBv, dBc);
  DIV_(eps, eps, x);
  for ( i = 0; i < npt; i++ ) {
    /* vc[i] *= eps; */
    MUL_X_(vc[i], eps);
    /* crl[i] += (fr[i] + 1) * vc[i] */
    ADD_SI_(x, fr[i], 1);
    FMA_X_(crl[i], x, vc[i]);
  }
  FMA_(B, dBv, eps, Bv0);
  CLEAR_(dBc);
  CLEAR_(dBv);
  CLEAR_(x);
}



/* Sum_{m = 3 to infinity} { rho^m [c(k)]^m }_{l+2} / m */
__inline static void get_ksum(mpfr_t s, int l, int npt,
    mpfr_t **ck, mpfr_t *w, mpfr_t s2)
{
  int i, j, k, m;
  mpfr_t *a, x, y, z;

  INIT_(x);
  INIT_(y);
  INIT_(z);
  if ( s2 != NULL )
    SET_SI_(s2, 0);
  SET_SI_(s, 0);
  MAKE1DARR(a, l + 1);
  for ( i = 0; i < npt; i++ ) {
    SET_SI_(y, 0);
    SET_SI_(z, 0);
    for ( j = 0; j <= l; j++ )
      SET_(a[j], ck[j][i]);
    for ( m = 2; m <= l + 2; m++ ) { /* rho^m [c(k)]^m */
      /* update a[l + 2 - m], ..., a[0] */
      for ( j = l + 2 - m; j >= 0; j-- ) {
        /* a'[j] = Sum_{k = 0 to j} a[k] ck[j - k]
         * we start from large j to preserve small-index data */
        /* a[j] = a[j] * ck[0][i]; */
        MUL_X_(a[j], ck[0][i]);
        for ( k = j - 1; k >= 0; k-- )
          /* a[j] += a[k] * ck[j - k][i]; */
          FMA_X_(a[j], a[k], ck[j - k][i]);
      }
      if ( m == 2 ) continue;
      /*
      y += a[l + 2 - m] / m;
      z += a[l + 2 - m] * (m - 2) / (2*m);
      */
      DIV_SI_(x, a[l + 2 - m], m);
      ADD_X_(y, x);
      DIV_SI_X_(x, 2);
      MUL_SI_X_(x, m - 2);
      ADD_X_(z, x);
    }
    FMA_X_(s, w[i], y);
    if ( s2 != NULL ) {
      FMA_X_(s2, w[i], z);
    }
  }
  FREE1DARR(a, l + 1);
  CLEAR_(x);
  CLEAR_(y);
  CLEAR_(z);
}



/* save the header for the virial file */
__inline static char *savevirhead(const char *fn, const char *title,
    int dim, int nmax, int dohnc, int mkcorr, int npt, double rmax,
    clock_t inittime)
{
  FILE *fp;
  int prec = (int) mpfr_get_default_prec();
  static char fndef[256];

  if ( fn == NULL ) {
    sprintf(fndef, "%sBn%s%sD%dn%dR%.0fM%dp%d.dat",
        title ? title : "", dohnc ? "HNC" : "PY", mkcorr ? "c" : "",
        dim, nmax, rmax, npt, prec);
    fn = fndef;
  }
  xfopen(fp, fn, "w", return NULL);
  fprintf(fp, "# %s %s %d %.14f %d %d | n Bc Bv Bm [Bh Br | corr] | %.3fs\n",
      dohnc ? "HNC" : "PY", mkcorr ? "corr" : "",
      nmax, rmax, npt, prec, (double) inittime / CLOCKS_PER_SEC);
  fclose(fp);
  return (char *) fn;
}



/* print a virial coefficient */
__inline static int printB(const char *name, int dim, int n,
    mpfr_t B, mpfr_t B2p, mpfr_t volp, const char *ending)
{
  mpfr_t x, xabs;
  int px = 0;

  if ( CMP_SI_(B, 0) == 0 ) return 0;
  mpfr_printf("%s(%3d) = %15.7Re", name, n, B);
  if ( CMP_SI_(B2p, 0) != 0 ) {
    INIT_(x);
    INIT_(xabs);
    DIV_(x, B, B2p);
    mpfr_printf(" (%15.8Re", x);
    DIV_(x, B, volp);
    ABS_(xabs, x);
    px = ( CMP_SI_(xabs, 10000) < 0 && dim % 2 == 1 );
    if ( !px ) mpfr_printf(")");
    else mpfr_printf(",%+12.6Rf)", x);
    CLEAR_(x);
    CLEAR_(xabs);
  }
  printf("%s", ending);
  return px;
}



/* save a virial coefficient to file */
__inline static void saveB(FILE *fp, mpfr_t B, mpfr_t B2p)
{
  if ( CMP_SI_(B, 0) == 0 ) {
    fprintf(fp, " 0");
  } else if ( CMP_SI_(B2p, 0) == 0 ) {
    mpfr_fprintf(fp, "%+24.14Re", B);
  } else {
    mpfr_t x;
    INIT_(x);
    DIV_(x, B, B2p);
    mpfr_fprintf(fp, "%+24.14Re", x);
    CLEAR_(x);
  }
}



/* save virial coefficients */
__inline static int savevir(const char *fn, int dim, int n,
    mpfr_t Bc, mpfr_t Bv, mpfr_t Bm, mpfr_t Bh, mpfr_t Br,
    mpfr_t By, mpfr_t B2, int mkcorr, mpfr_t fcorr)
{
  FILE *fp;
  mpfr_t B2p, B2y, volp, voly, x;

  INIT_(B2p);
  INIT_(B2y);
  INIT_(volp);
  INIT_(voly);
  INIT_(x);
  POW_SI_(B2p, B2, n - 1);
  POW_SI_(B2y, B2, n - 2);
  /* print the result on screen */
  /* volp = B2p / pow(2, (n-1)*(dim-1)); */
  SET_SI_(x, 2);
  POW_SI_(volp, x, (n-1)*(dim-1));
  DIV_(volp, B2p, volp);
  POW_SI_(voly, x, (n-2)*(dim-1));
  DIV_(voly, B2y, voly);
  printB("Bc", dim, n, Bc, B2p, volp, ", ");
  printB("Bv", dim, n, Bv, B2p, volp, ", ");
  /* when making corrections, Bm is the corrected value */
  printB("Bm", dim, n, Bm, B2p, volp, "");
  if ( mkcorr ) {
    mpfr_printf(", %9.6Rf", fcorr);
  } else { /* the following are useless when making corrections */
    if ( CMP_SI_(Bm, 0) == 0 ) printf("\n");
    printB("Bh", dim, n, Bh, B2p, volp, ", ");
    printB("Br", dim, n, Br, B2p, volp, "");
    printB("By", dim, n - 1, By, B2y, voly, ", ");
  }
  printf("\n");

  if (fn != NULL) {
    xfopen(fp, fn, "a", return -1);
    fprintf(fp, "%4d", n);
    /* normalize the virial coefficients */
    if ( mkcorr ) {
      saveB(fp, Bc, B2p);
      saveB(fp, Bv, B2p);
      saveB(fp, Bm, B2p);
      mpfr_fprintf(fp, " %+20.14Rf", fcorr);
    } else {
      saveB(fp, Bc, B2p);
      saveB(fp, Bv, B2p);
      saveB(fp, Bm, B2p);
      saveB(fp, Bh, B2p);
      saveB(fp, Br, B2p);
    }
    fprintf(fp, " |");
    saveB(fp, By, B2y);
    fprintf(fp, "\n");
    fclose(fp);
  }
  CLEAR_(B2p);
  CLEAR_(B2y);
  CLEAR_(volp);
  CLEAR_(voly);
  CLEAR_(x);

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
  fprintf(stderr, "virial coefficients saved to %s\n", fn);
  return 0;
}



/* save c(r) or t(r) file */
__inline static int savecrtr(const char *fn, int l, int npt,
    mpfr_t *ri, mpfr_t *cr, mpfr_t *tr, mpfr_t *vc, mpfr_t *yrl)
{
  FILE *fp;
  int i;

  if (fn == NULL) return -1;
  xfopen(fp, fn, (l == 1) ? "w" : "a", return -1);
  for ( i = 0; i < npt; i++ ) {
    mpfr_fprintf(fp, "%10.7Rf %20.12Rf %20.12Rf %d",
        ri[i], cr[i], tr[i], l);
    if ( vc != NULL )
      mpfr_fprintf(fp, " %20.12Rf", vc[i]);
    if ( yrl != NULL )
      mpfr_fprintf(fp, " %20.12Rf", yrl[i]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  fclose(fp);
  return 0;
}



/* compute mu = Int (c - h t / 2) dr */
__inline static void get_Bm_singer(mpfr_t s, int l, int npt,
    mpfr_t **cr, mpfr_t **tr, mpfr_t *w)
{
  int i, u;
  mpfr_t x, y;

  INIT_(x);
  INIT_(y);
  for ( i = 0; i < npt; i++ ) {
    SET_(x, cr[l][i]);
    for ( u = 0; u < l; u++ ) {
      /* x -= tr[l-u][i] * (cr[u][i] + tr[u][i]) / 2; */
      FAM_(y, cr[u][i], tr[u][i], tr[l-u][i]);
      DIV_SI_X_(y, 2);
      SUB_X_(x, y);
    }
    FMA_X_(s, x, w[i]);
  }
  CLEAR_(x);
  CLEAR_(y);
  MUL_SI_X_(s, -(l+1));
  DIV_SI_X_(s, (l+2));
}




/* compute -beta F = (1/2) Int (c - c t - t^2 / 2) dr (real part) */
__inline static void get_Bh_singer(mpfr_t s, int l, int npt,
    mpfr_t **cr, mpfr_t **tr, mpfr_t *w)
{
  int i, u;
  mpfr_t x, y;

  INIT_(x);
  INIT_(y);
  /* compute -2 beta F */
  SET_SI_(s, 0);
  for ( i = 0; i < npt; i++ ) {
    SET_(x, cr[l][i]);
    for ( u = 0; u < l; u++ ) {
      /* x -= tr[l-u][i] * (cr[u][i] + tr[u][i] / 2); */
      DIV_SI_(y, tr[u][i], 2);
      FAM_(y, y, cr[u][i], tr[l-u][i]);
      SUB_X_(x, y);
    }
    FMA_X_(s, x, w[i]);
  }
  CLEAR_(x);
  CLEAR_(y);
  MUL_SI_X_(s, -(l+1));
  DIV_SI_X_(s, 2);
}



/* d^2_rho beta P = Int (c + h t) dr */
__inline static void get_Bx_py(mpfr_t s, int l, int npt,
    mpfr_t **cr, mpfr_t **tr, mpfr_t *w)
{
  int i, u;
  mpfr_t x, y;

  INIT_(x);
  INIT_(y);
  SET_SI_(s, 0);
  for ( i = 0; i < npt; i++ ) {
    SET_(x, cr[l][i]);
    for ( u = 0; u < l; u++ ) {
      ADD_(y, cr[u][i], tr[u][i]);
      FMA_X_(x, tr[l-u][i], y);
    }
    FMA_X_(s, x, w[i]);
  }
  CLEAR_(x);
  CLEAR_(y);
  DIV_SI_X_(s, -(l+1)*(l+2));
}



/* d^2_rho beta P = (-1/2) Int (c - t c) dr (real part) */
__inline static void get_Bp_py(mpfr_t s, int l, int npt,
    mpfr_t **cr, mpfr_t **tr, mpfr_t *w)
{
  int i, u;
  mpfr_t x, y;

  INIT_(x);
  INIT_(y);
  //NEG_(s, s);
  SET_SI_(s, 0);
  for ( i = 0; i < npt; i++ ) {
    SET_(x, cr[l][i]);
    for ( u = 0; u < l; u++ ) {
      /* x -= tr[l-u][i] * cr[u][i]; */
      MUL_(y, tr[l-u][i], cr[u][i]);
      SUB_X_(x, y);
    }
    FMA_X_(s, x, w[i]);
  }
  CLEAR_(x);
  CLEAR_(y);
  DIV_SI_X_(s, -2);
}



#endif

