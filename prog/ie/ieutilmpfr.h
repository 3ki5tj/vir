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



#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (int) (n); i_++) { \
    INIT_(arr[i_]); \
    SET_SI_(arr[i_], 0); } }

#define FREE1DARR(arr, n) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) \
    CLEAR_(arr[i_]); \
    free(arr); }

#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr, n1, n2) { \
  FREE1DARR(arr[0], (n1) * (n2)); \
  free(arr); }

#define COPY1DARR(a, b, n) { int i_; \
  for (i_ = 0; i_ < (int) (n); i_++) \
    SET_(a[i_], b[i_]); }



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



/* compute t(k) from c(k) from the Ornstein-Zernike relation */
__inline static void get_tk_oz(int l, int npt, mpfr_t **ck, mpfr_t **tk)
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
static void get_yr_hnc(int l, int nmax, int npt,
    mpfr_t **yr, mpfr_t *trl)
{
  int i, j, k, jl;
  mpfr_t x, powtr, *yro;

  INIT_(x);
  INIT_(powtr);
  MAKE1DARR(yro, nmax - 1);
  /* y(r) = exp(t(r))
   * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
  /* yr = yr0 * exp(rho^l t_l)
   *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
   *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
   * y[l...n](r) are updated */
  for ( i = 0; i < npt; i++) {
    for ( j = 0; j < nmax - 1; j++ )
      SET_(yro[j], yr[j][i]);
    SET_SI_(powtr, 1);
    for ( j = 1; j*l < nmax - 1; j++ ) {
      /* powtr = tl^j/j!; powtr *= tl/j; */
      DIV_SI_(x, trl[i], j);
      MUL_X_(powtr, x);
      for ( jl = j * l, k = 0; k + jl < nmax - 1; k++ ) {
        /* yr_{k + jl} +=  yr0_k * tl^j/j! */
        FMA_X_(yr[jl+k][i], yro[k], powtr);
      }
    }
  }
  CLEAR_(x);
  CLEAR_(powtr);
  FREE1DARR(yro, nmax - 1);
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
    SET_(B, Bc0);
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
    int dim, int nmax, int doHNC, int mkcorr, int npt, double rmax)
{
  FILE *fp;
  static char fndef[256];

  if ( fn == NULL ) {
    sprintf(fndef, "%sBn%s%sD%dn%dR%.0fM%d.dat",
        title ? title : "", doHNC ? "HNC" : "PY", mkcorr ? "c" : "",
        dim, nmax, rmax, npt);
    fn = fndef;
  }
  xfopen(fp, fn, "w", return NULL);
  fprintf(fp, "# %s %s %d %.14f %d | n Bc Bv Bm [Bh Br | corr]\n",
      doHNC ? "HNC" : "PY", mkcorr ? "corr" : "",
      nmax, rmax, npt);
  fclose(fp);
  return (char *) fn;
}



/* print a virial coefficient */
__inline static void printB(const char *name, int n, mpfr_t B,
    mpfr_t B2p, mpfr_t volp, const char *ending)
{
  if ( B == 0 ) return;
  mpfr_printf("%s(%3d) = %15.8Re", name, n, B);
  if ( B2p != 0 ) {
    mpfr_t x, y;
    INIT_(x);
    INIT_(y);
    DIV_(x, B, B2p);
    DIV_(y, B, volp);
    mpfr_printf(" (%15.8Re, %11.6Rf)", x, y);
    CLEAR_(x);
    CLEAR_(y);
  }
  printf("%s", ending);
}



/* save virial coefficients */
__inline static int savevir(const char *fn, int dim, int l,
    mpfr_t Bc, mpfr_t Bv, mpfr_t Bm, mpfr_t Bh, mpfr_t Br,
    mpfr_t B2p, int mkcorr, mpfr_t fcorr)
{
  FILE *fp;
  mpfr_t volp, x;

  INIT_(volp);
  INIT_(x);
  /* print the result on screen */
  /* volp = B2p / pow(2, (l+1)*(dim-1)); */
  SET_SI_(x, 2);
  POW_SI_(x, x, (l+1)*(dim-1));
  DIV_(volp, B2p, x);
  printB("Bc", l+2, Bc, B2p, volp, ", ");
  printB("Bv", l+2, Bv, B2p, volp, ", ");
  /* when making corrections, Bm is the corrected value */
  printB("Bm", l+2, Bm, B2p, volp, "");
  if ( mkcorr ) {
    mpfr_printf(", %9.6Rf\n", fcorr);
  } else { /* the following are useless when making corrections */
    printf("\n");
    if ( !mpfr_zero_p(Bh) || !mpfr_zero_p(Br) ) {
      printB("Bh", l+2, Bh, B2p, volp, ", ");
      printB("Br", l+2, Br, B2p, volp, "\n");
    }
  }

  if (fn != NULL) {
    mpfr_t Xc, Xv, Xm, Xh, Xr;
    xfopen(fp, fn, "a", return -1);
    INIT_(Xc);
    INIT_(Xv);
    INIT_(Xm);
    INIT_(Xh);
    INIT_(Xr);
    /* normalize the virial coefficients */
    if ( B2p != 0 ) {
      DIV_(Xc, Bc, B2p);
      DIV_(Xv, Bv, B2p);
      DIV_(Xm, Bm, B2p);
      DIV_(Xh, Bh, B2p);
      DIV_(Xr, Br, B2p);
    } else {
      SET_(Xc, Bc);
      SET_(Xv, Bv);
      SET_(Xm, Bm);
      SET_(Xh, Bh);
      SET_(Xr, Br);
    }
    if ( mkcorr ) {
      mpfr_fprintf(fp, "%4d%+24.14Re%+24.14Re%+24.14Re %+18.14Rf\n",
          l + 2, Bc, Bv, Bm, fcorr);
    } else {
      mpfr_fprintf(fp, "%4d%+24.14Re%24.14Re%24.14Re%24.14Re%24.14Re\n",
          l + 2, Bc, Bv, Bm, Bh, Br);
    }
    fclose(fp);
    CLEAR_(Xc);
    CLEAR_(Xv);
    CLEAR_(Xm);
    CLEAR_(Xh);
    CLEAR_(Xr);
  }
  CLEAR_(volp);
  CLEAR_(x);

  return 0;
}



/* save c(r) or t(r) file */
__inline static int savecrtr(const char *fn, int l, int npt,
    mpfr_t *ri, mpfr_t *cr, mpfr_t *tr, mpfr_t *vc, mpfr_t **yr)
{
  FILE *fp;
  int i, j;

  if (fn == NULL) return -1;
  xfopen(fp, fn, (l == 1) ? "w" : "a", return -1);
  for ( i = 0; i < npt; i++ ) {
    mpfr_fprintf(fp, "%10.7Rf %20.12Rf %20.12Rf %d",
        ri[i], cr[i], tr[i], l);
    if ( vc != NULL )
      mpfr_fprintf(fp, " %20.12Rf", vc[i]);
    if ( yr != NULL )
      for ( j = 1; j <= l; j++ )
        mpfr_fprintf(fp, " %20.12Rf", yr[j][i]);
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

