/* Computing the virial coefficients of a odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ieodmp.c -lmpfr
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpfft.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#ifdef D
int dim = D;
int K = (D - 1)/2;
#else
int dim = 3;
int K = 1;
#endif

int prec = 0;
int nmax = 10;
char *dr = NULL;
int numpt = 32768;
int doHNC = 0;
int mkcorr = 0;



#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (n); i_++) { \
    mpfr_init(arr[i_]); \
    mpfr_set_si(arr[i_], 0, MPFR_RNDN); } }

#define FREE1DARR(arr, n) { int i_; \
  for (i_ = 0; i_ < (n); i_++) \
    mpfr_clear(arr[i_]); \
  free(arr); }

#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); \
  MAKE1DARR(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) \
    arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr, n1, n2) { \
  FREE1DARR(arr[0], (n1) * (n2)); \
  free(arr); }



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the PY/HNC closure for odd-dimensional hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-b", "%s", &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-p", "%d", &prec, "float-point precision in bits");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if (dr == NULL) dr = "0.001";
  if (dim < 3 || dim % 2 == 0) argopt_help(ao);
  K = (dim - 1)/2;
  if (prec == 0) prec = 128 * dim;
  argopt_dump(ao);
  argopt_close(ao);
}



/* get the coefficient c_l of spherical Bessel function jn(x)
 *  jn(x) = Sum_{l = 0 to n} c_l [sin(x) or cos(x)] / x^{n + l + 1},
 * where the l'th term is sin(x) if n + l is even, or cos(x) otherwise */
static void getjn(int *c, int n)
{
  int i, k;
  const char *fs[2] = {"sin(x)", "cos(x)"};

  c[0] = 1; /* j0 = sin(x)/x; */
  /* j_n(x)/x^n = (-1/x d/dx)^n j_0(x) */
  for ( k = 1; k <= n; k++ ) { /* k'th round */
    c[k] = 0;
    for ( i = k - 1; i >= 0; i-- ) {
      c[i + 1] += c[i] * (k + i);
      c[i] *= 1 - (k + i) % 2 * 2;
    }
  }
  printf("j%d(x) =", n);
  for ( i = 0; i <= n; i++ )
    printf(" %+d*%s/x^%d", c[i], fs[(i+n)%2], n + i + 1);
  printf("\n");
}



/* compute the D-dimensional Fourier transform of a spherical transform
 *    out(k) = 2 fac0 Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' is the intermediate array
 * */
static void sphr11(int npt, mpfr_t *in, mpfr_t *out,
    mpfr_t dx, mpfr_t fac0, mpfr_t *arr,
    int *coef, mpfr_t **r2p, mpfr_t **k2q)
{
  int i, l, iscos;
  mpfr_t fac, x;

  mpfr_init(fac);
  mpfr_init(x);

  /* fac = fac0 * dx */
  mpfr_mul(fac, fac0, dx, MPFR_RNDN);

  /* clear the output */
  for ( i = 0; i < npt; i++ )
    mpfr_set_si(out[i], 0, MPFR_RNDN); /* out[i] = 0; */

  /* several rounds of transforms */
  for ( l = 0; l < K; l++ ) {
    /* decide if we want to do a sine or cosine transform */
    iscos = (K + l + 1) % 2;
    for ( i = 0; i < npt; i++ ) {
      /*
       * arr[i] = in[i] * pow(dx*(2*i + 1)/2, K - l) * coef[l];
       * or
       * arr[i] = in[i] * coef[l] * r2p[l][i];
       * */
      mpfr_mul_si(x, r2p[l][i], coef[l], MPFR_RNDN);
      mpfr_mul(arr[i], in[i], x, MPFR_RNDN);
    }

    /* do the sine or cosine transform */
    if ( iscos ) {
      mpcost11(arr, npt, NULL);
    } else {
      mpsint11(arr, npt, NULL);
    }

    for ( i = 0; i < npt; i++ ) {
      /*
       * out[i] += arr[i] * fac / pow(dk*(2*i + 1)/2, K + l);
       * or
       * out[i] += arr[i] * fac * k2q[l][i];
       * */
      mpfr_mul(x, fac, k2q[l][i], MPFR_RNDN);
      mpfr_mul(x, x, arr[i], MPFR_RNDN);
      mpfr_add(out[i], out[i], x, MPFR_RNDN);
    }
  }

  mpfr_clear(fac);
  mpfr_clear(x);
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, const char *sdr, int doHNC)
{
  mpfr_t dr, dk, B2, pi2, facr2k, fack2r, tmp1, tmp2;
  mpfr_t Bc0, dBc, Bv0, dBv, eps;
  mpfr_t **ck, **tk, *fk, *fr, *Bc, *Bv;
  mpfr_t *arr, *tlk, *tlr, *clk, *clr;
  mpfr_t *powtlr, **yr0, **yr, *vc;
  mpfr_t **r2p, **invr2p, **k2p, **invk2p, *r2Dm1;
  int i, dm, l, u, j, jl, k;
  int *coef;

  mpfr_init(dr);
  mpfr_init(dk);
  mpfr_init(B2);
  mpfr_init(pi2);
  mpfr_init(facr2k);
  mpfr_init(fack2r);
  mpfr_init(tmp1);
  mpfr_init(tmp2);

  mpfr_init(Bc0);
  mpfr_init(dBc);
  mpfr_init(Bv0);
  mpfr_init(dBv);
  mpfr_init(eps);

  mpfr_set_str(dr, sdr, 10, MPFR_RNDN);
  /* dk = PI/npt/dr */
  mpfr_const_pi(dk, MPFR_RNDN);
  mpfr_div_si(dk, dk, npt, MPFR_RNDN);
  mpfr_div(dk, dk, dr, MPFR_RNDN);

  /* the bin index at which the hard core stops */
  dm = (int) (1/atof(sdr) + .5);

  /* print out the basic information */
  {
    double rmax = atof(sdr) * npt;
    double dk_d = mpfr_get_d(dk, MPFR_RNDN);

    i = (int) mpfr_get_default_prec();
    printf("precision %d, rmax %g, dk %g, %d bins in the hard core\n",
        i, rmax, dk_d, dm);
  }

  /* compute the coefficients of the spherical Bessel function */
  xnew(coef, K);
  getjn(coef, K - 1);

  MAKE1DARR(arr, npt);

  MAKE2DARR(r2p, K, npt) /* r^{K - l} for l = 0, ..., K */
  MAKE2DARR(k2p, K, npt) /* k^{K - l} for l = 0, ..., K */
  MAKE2DARR(invr2p, K, npt) /* r^{-K-l} for l = 0, ..., K */
  MAKE2DARR(invk2p, K, npt) /* k^{-K-l} for l = 0, ..., K */
  {
    mpfr_t ri, rl, invrl, ki, kl, invkl;

    mpfr_init(ri);
    mpfr_init(rl);
    mpfr_init(invrl);
    mpfr_init(ki);
    mpfr_init(kl);
    mpfr_init(invkl);
    for ( i = 0; i < npt; i++ ) {
      /* ri = dr * (i*2 + 1)/2; */
      mpfr_mul_si(ri, dr, i*2 + 1, MPFR_RNDN);
      mpfr_div_si(ri, ri, 2, MPFR_RNDN);
      mpfr_set_si(rl, 1, MPFR_RNDN); /* rl = 1 */
      /* ki = dk * (i*2 + 1)/2; */
      mpfr_mul_si(ki, dk, i*2 + 1, MPFR_RNDN);
      mpfr_div_si(ki, ki, 2, MPFR_RNDN);
      mpfr_set_si(kl, 1, MPFR_RNDN); /* kl = 1 */
      for ( l = 1; l <= K; l++ ) {
        mpfr_mul(rl, rl, ri, MPFR_RNDN); /* rl *= ri */
        mpfr_set(r2p[K-l][i], rl, MPFR_RNDN); /* r2p[K-l][i] = rl; */
        mpfr_mul(kl, kl, ki, MPFR_RNDN); /* kl *= ki */
        mpfr_set(k2p[K-l][i], kl, MPFR_RNDN); /* k2p[K-l][i] = kl; */
      }
      mpfr_si_div(invrl, 1, rl, MPFR_RNDN); /* invrl = 1/rl; */
      mpfr_si_div(invkl, 1, kl, MPFR_RNDN); /* invkl = 1/kl; */
      for ( l = 0; l < K; l++ ) {
        mpfr_set(invr2p[l][i], invrl, MPFR_RNDN); /* invr2p[l][i] = invrl; */
        mpfr_set(invk2p[l][i], invkl, MPFR_RNDN); /* invk2p[l][i] = invkl; */
        mpfr_div(invrl, invrl, ri, MPFR_RNDN); /* invrl /= ri; */
        mpfr_div(invkl, invkl, ki, MPFR_RNDN); /* invkl /= ki; */
      }
    }
    mpfr_clear(ri);
    mpfr_clear(rl);
    mpfr_clear(invrl);
    mpfr_clear(ki);
    mpfr_clear(kl);
    mpfr_clear(invkl);
  }

  /* compute r^(D - 1) dr */
  MAKE1DARR(r2Dm1, npt);
  for (i = 0; i < npt; i++) {
    /* tmp1 = dr * (i*i + 1)/2; */
    mpfr_mul_si(tmp1, dr, i*2+1, MPFR_RNDN);
    mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
    /* tmp2 = tmp1^(dim - 1) * dr */
    mpfr_pow_si(tmp2, tmp1, dim - 1, MPFR_RNDN);
    /* r2Dm1[i] = pow(tmp1, dim - 1) * dr; */
    mpfr_mul(r2Dm1[i], tmp2, dr, MPFR_RNDN);
  }

  /* B2 = (PI*2)^K/(2 K + 1)!!
   * facr2k = (PI*2)^K */
  mpfr_set_si(B2, 1, MPFR_RNDN); /* B2 = 1 */
  mpfr_set_si(facr2k, 1, MPFR_RNDN); /* facr2k = 1 */
  mpfr_const_pi(pi2, MPFR_RNDN);
  mpfr_mul_si(pi2, pi2, 2, MPFR_RNDN); /* pi2 = PI*2 */
  for (i = 1; i <= K; i++) {
    /* B2 *= PI*2/(i*2 + 1); */
    mpfr_mul(B2, B2, pi2, MPFR_RNDN);
    mpfr_div_si(B2, B2, i*2 + 1, MPFR_RNDN);
    /* facr2k *= PI*2; */
    mpfr_mul(facr2k, facr2k, pi2, MPFR_RNDN);
  }
  /* fack2r = 1/(facr2k * PI*2); */
  mpfr_si_div(fack2r, 1, facr2k, MPFR_RNDN);
  mpfr_div(fack2r, fack2r, pi2, MPFR_RNDN);

  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  /* Bc[2] = Bv[2] = B2; */
  mpfr_set(Bc[2], B2, MPFR_RNDN);
  mpfr_set(Bv[2], B2, MPFR_RNDN);

  /* construct f(r) and f(k) */
  MAKE1DARR(fr, npt);
  MAKE1DARR(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    mpfr_set_si(fr[i], (i < dm) ? -1 : 0, MPFR_RNDN);
  /* f(r) --> f(k) */
  sphr11(npt, fr, fk, dr, facr2k, arr, coef, r2p, invk2p);

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ )
    mpfr_set(ck[0][i], fk[i], MPFR_RNDN); /* ck[0][i] = fk[i]; */

  MAKE1DARR(clk, npt); MAKE1DARR(clr, npt);
  MAKE1DARR(tlk, npt); MAKE1DARR(tlr, npt);

  if ( doHNC ) {
    MAKE1DARR(powtlr, npt);
    MAKE2DARR(yr0, nmax + 1, npt);
    MAKE2DARR(yr, nmax + 1, npt);
    for ( i = 0; i < npt; i++ ) {
      mpfr_set_si(yr0[0][i], 1, MPFR_RNDN);
      mpfr_set_si(yr[0][i], 1, MPFR_RNDN);
    }
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute t_l(k) */
    for ( i = 0; i < npt; i++ )
      mpfr_set_si(tlk[i], 0, MPFR_RNDN);
    for ( u = 0; u < l; u++ )
      for ( i = 0; i < npt; i++ ) {
        /* tlk[i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]); */
        mpfr_add(tmp1, ck[u][i], tk[u][i], MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, ck[l-1-u][i], MPFR_RNDN);
        mpfr_add(tlk[i], tlk[i], tmp1, MPFR_RNDN);
      }
    for ( i = 0; i < npt; i++ )
      mpfr_set(tk[l][i], tlk[i], MPFR_RNDN); /* tk[l][i] = tlk[i]; */

    /* t_l(k) --> t_l(r) */
    sphr11(npt, tlk, tlr, dk, fack2r, arr, coef, k2p, invr2p);

    if ( doHNC ) {
      /* Hypernetted chain approximation
       * y(r) = exp(t(r))
       * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
      /* back up, only need to save yr[1..nmax-l] */
      for ( u = 1; u <= nmax - l; u++ )
        for ( i = 0; i < npt; i++ )
          mpfr_set(yr0[u][i], yr[u][i], MPFR_RNDN);
      /* yr = yr0 * exp(rho^l t_l)
       *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
       *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
       * y[l...n](r) are updated */
      for ( i = 0; i < npt; i++ )
        mpfr_set_si(powtlr[i], 1, MPFR_RNDN); /* powtlr[i] = 1; */
      for ( j = 1; j*l <= nmax; j++ ) {
        /* powtl(r) = tl(r)^j/j! */
        for ( i = 0; i < npt; i++) {
          mpfr_div_si(tmp1, tlr[i], j, MPFR_RNDN);
          mpfr_mul(powtlr[i], powtlr[i], tmp1, MPFR_RNDN);
          /* powtlr[i] *= tlr[i]/j; */
        }
        for ( jl = j * l, k = 0; k + jl <= nmax; k++ )
          /* yr_{k + jl} +=  yr0_k * tl^j/j! */
          for ( i = 0; i < npt; i++ ) {
            mpfr_mul(tmp1, yr0[k][i], powtlr[i], MPFR_RNDN);
            mpfr_add(yr[jl+k][i], yr[jl+k][i], tmp1, MPFR_RNDN);
            /* yr[jl+k][i] += yr0[k][i] * powtlr[i]; */
          }
      }

      /* compute c(r) = (1 + f) y(r) - (1 + t(r)) */
      for ( i = 0; i < dm; i++ )
        mpfr_neg(clr[i], tlr[i], MPFR_RNDN); /* clr[i] = -tlr[i]; */
      for ( i = dm; i < npt; i++ ) {
        mpfr_sub(clr[i], yr[l][i], tlr[i], MPFR_RNDN);
        /* clr[i] = yr[l][i] - tlr[i]; */
      }

      if ( mkcorr && l >= 2 ) { /* make a correction */
        /* `vc' is the trial correction function to y(r) */
        for ( i = 0; i < npt; i++ ) {
          /* For the hard-sphere fluid, the trial function cl(r)
           *   is equivalent to yl(r) - tl(r)
           * By definition, we have c(r) = [1 + f(r)] y(r) - [1 + t(r)],
           *   and f(r) = (r <= 1) ? -1 : 0
           * So
           *   c_l(r) = y_l(r) - t_l(r),  r > 1
           *          = -t_l(r),          r <= 1
           * Note that the r <= 1 branch is cancelled by 1 + f(r) */
          mpfr_set(vc[i], clr[i], MPFR_RNDN);
        }

        mpfr_set_si(Bc0, 0, MPFR_RNDN);
        mpfr_set_si(dBc, 0, MPFR_RNDN);
        for ( i = 0; i < npt; i++ ) {
          /* Bc0 += clr[i] * r2Dm1[i]; */
          mpfr_mul(tmp1, clr[i], r2Dm1[i], MPFR_RNDN);
          mpfr_add(Bc0, Bc0, tmp1, MPFR_RNDN);

          /* dBc += vc[i] * (1 + fr[i]) * r2Dm1[i]; */
          mpfr_add_si(tmp1, fr[i], 1, MPFR_RNDN);
          mpfr_mul(tmp1, tmp1, r2Dm1[i], MPFR_RNDN);
          mpfr_mul(tmp1, tmp1, vc[i], MPFR_RNDN);
          mpfr_add(dBc, dBc, tmp1, MPFR_RNDN);
        }
        /* tmp1 = -2*B2*dim/(l+2); */
        mpfr_mul_si(tmp1, B2, -dim*2, MPFR_RNDN);
        mpfr_div_si(tmp1, tmp1, l+2, MPFR_RNDN);
        mpfr_mul(Bc0, Bc0, tmp1, MPFR_RNDN); /* Bc0 *= tmp1; */
        mpfr_mul(dBc, dBc, tmp1, MPFR_RNDN); /* dBc *= tmp1; */

        /* Bv0 = B2 * (yr[l][dm] + yr[l][dm-1])/2; */
        mpfr_add(tmp1, yr[l][dm], yr[l][dm-1], MPFR_RNDN);
        mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
        mpfr_mul(Bv0, B2, tmp1, MPFR_RNDN);
        /* dBv = B2 * (vc[dm] + vc[dm-1])/2; */
        mpfr_add(tmp1, vc[dm], vc[dm-1], MPFR_RNDN);
        mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
        mpfr_mul(dBv, B2, tmp1, MPFR_RNDN);
        /* eps = -(Bv0 - Bc0) / (dBv - dBc); */
        mpfr_sub(tmp1, Bc0, Bv0, MPFR_RNDN);
        mpfr_sub(tmp2, dBv, dBc, MPFR_RNDN);
        mpfr_div(eps, tmp1, tmp2, MPFR_RNDN);
        for ( i = 0; i < npt; i++ ) {
          mpfr_mul(tmp1, eps, vc[i], MPFR_RNDN);
          mpfr_add(yr[l][i], yr[l][i], tmp1, MPFR_RNDN);
          /* yr[l][i] += eps * vc[i]; */
          mpfr_add_si(tmp2, fr[i], 1, MPFR_RNDN);
          mpfr_mul(tmp2, tmp2, tmp1, MPFR_RNDN);
          mpfr_add(clr[i], clr[i], tmp2, MPFR_RNDN);
          /* clr[i] += (fr[i] + 1) * eps * vc[i] */
        }
      }
      /* Bv[l+2] = B2*(yr[l][dm] + yr[l][dm-1])/2; */
      mpfr_add(tmp1, yr[l][dm], yr[l][dm-1], MPFR_RNDN);
      mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ )
        mpfr_neg(clr[i], tlr[i], MPFR_RNDN); /* clr[i] = -tlr[i]; */
      for ( i = dm; i < npt; i++ )
        mpfr_set_si(clr[i], 0, MPFR_RNDN); /* clr[i] = 0; */
      /* Bv[l+2] = B2*(tlr[dm] + tlr[dm-1])/2; */
      mpfr_add(tmp1, tlr[dm], tlr[dm-1], MPFR_RNDN);
      mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 2 (2 pi)^K/(2*K - 1)!! r^(2*K) dr
     *           = -2*B2*D/(l+2) Int c_l(r) r^(2*K) dr */
    mpfr_set_si(Bc[l+2], 0, MPFR_RNDN); /* Bc[l+2] = 0; */
    for ( i = 0; i < npt; i++ ) {
      /* Bc[l+2] += clr[i] * r2Dm1[i]; */
      mpfr_mul(tmp1, clr[i], r2Dm1[i], MPFR_RNDN);
      mpfr_add(Bc[l+2], Bc[l+2], tmp1, MPFR_RNDN);
    }
    /* Bc[l+2] *= -B2*dim*2/(l+2); */
    mpfr_mul_si(tmp1, B2, -dim*2, MPFR_RNDN);
    mpfr_div_si(tmp1, tmp1, l+2, MPFR_RNDN);
    mpfr_mul(Bc[l+2], Bc[l+2], tmp1, MPFR_RNDN);

    /* tmp1 = Bc/B2^(l+1), tmp2 = Bv/B2^(l+1) */
    mpfr_pow_si(tmp1, B2, l+1, MPFR_RNDN);
    mpfr_div(tmp2, Bv[l+2], tmp1, MPFR_RNDN);
    mpfr_div(tmp1, Bc[l+2], tmp1, MPFR_RNDN);
    mpfr_printf("Bc(%3d) = %20.10Re (%20.12Re), "
                "Bv(%3d) = %20.10Re (%20.12Re)\n",
                l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c_l(r) --> c_l(k) */
    sphr11(npt, clr, clk, dr, facr2k, arr, coef, r2p, invk2p);

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ )
      mpfr_set(ck[l][i], clk[i], MPFR_RNDN); /* ck[l][i] = clk[i]; */
  }
  FREE1DARR(arr, npt);
  FREE1DARR(clr, npt); FREE1DARR(clk, npt);
  FREE1DARR(tlr, npt); FREE1DARR(tlk, npt);
  FREE1DARR(fr, npt); FREE1DARR(fk, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(ck, nmax - 1, npt);
  if ( doHNC ) {
    FREE1DARR(powtlr, npt);
    FREE2DARR(yr0, nmax + 1, npt);
    FREE2DARR(yr, nmax + 1, npt);
    if ( mkcorr ) {
      FREE1DARR(vc, npt);
    }
  }
  FREE1DARR(Bc, nmax + 1); FREE1DARR(Bv, nmax + 1);
  free(coef);
  FREE1DARR(r2Dm1, npt);
  FREE2DARR(r2p, K, npt); FREE2DARR(invr2p, K, npt);
  FREE2DARR(k2p, K, npt); FREE2DARR(invk2p, K, npt);

  mpfr_clear(dr);
  mpfr_clear(dk);
  mpfr_clear(B2);
  mpfr_clear(pi2);
  mpfr_clear(facr2k);
  mpfr_clear(fack2r);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);

  mpfr_clear(Bc0);
  mpfr_clear(dBc);
  mpfr_clear(Bv0);
  mpfr_clear(dBv);
  mpfr_clear(eps);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mpfr_set_default_prec(prec);

  intgeq(nmax, numpt, dr, doHNC);

  MPFFT_ARR1D_FREE();
  mpfr_free_cache();
  ssdelall();
  return 0;
}
