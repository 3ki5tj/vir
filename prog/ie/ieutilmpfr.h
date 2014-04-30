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
__inline static void get_corr1_hnc_hs(int l, int npt, int dm,
    mpfr_t *ylr, mpfr_t *crl, mpfr_t *fr, mpfr_t *rDm1,
    mpfr_t B2, mpfr_t *vc)
{
  int i;
  mpfr_t Bc0, dBc, Bv0, dBv, eps, x;

  if ( l < 2 ) return;

  INIT_(Bc0);
  INIT_(dBc);
  INIT_(Bv0);
  INIT_(dBv);
  INIT_(eps);
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

  CLEAR_(Bc0);
  CLEAR_(dBc);
  CLEAR_(Bv0);
  CLEAR_(dBv);
  CLEAR_(eps);
  CLEAR_(x);
}



#endif

