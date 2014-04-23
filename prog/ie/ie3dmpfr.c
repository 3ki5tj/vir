/* Computing the virial coefficients of the 3D hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ie3dmp.c -lmpfr
 * */
#include <stdio.h>
#include <math.h>
#include "mpfft.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"


int prec = 512;
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

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-b", "%s", &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-p", "%d", &prec, "float-point precision in bits");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if (dr == NULL) dr = "0.001";
  argopt_dump(ao);
  argopt_close(ao);
}



#if 0
/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr00(int npt, mpfr_t *in, mpfr_t *out,
    mpfr_t dx, mpfr_t dk, mpfr_t fac0, mpfr_t *arr)
{
  int i;
  mpfr_t fac, x;

  mpfr_init(fac);
  mpfr_init(x);
  for ( i = 1; i < npt; i++ ) { /* form in(x) * x */
    mpfr_mul_si(x, dx, i, MPFR_RNDN);
    mpfr_mul(arr[i], in[i], x, MPFR_RNDN);
    /* arr[i] = in[i] * x; */
  }

  /* do the sine transform */
  mpsint00(arr, npt, NULL);

  /* 1. the following dx is from integration
   * 2. divided by dk is for 1/k */
  /* fac = fac0 * dx / dk */
  mpfr_mul(fac, fac0, dx, MPFR_RNDN);
  mpfr_div(fac, fac, dk, MPFR_RNDN);

  for ( i = 1; i < npt; i++ ) { /* form out(k) / k */
    mpfr_div_si(x, fac, i, MPFR_RNDN);
    mpfr_mul(out[i], arr[i], x, MPFR_RNDN);
    /* out[i] = arr[i] * fac / i */
  }
  mpfr_clear(fac);
  mpfr_clear(x);
}
#endif



/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr11(int npt, mpfr_t *in, mpfr_t *out,
    mpfr_t dx, mpfr_t dk, mpfr_t fac0, mpfr_t *arr)
{
  int i;
  mpfr_t fac, x;

  mpfr_init(fac);
  mpfr_init(x);
  for ( i = 0; i < npt; i++ ) { /* form in(x) * x */
    mpfr_set_si(x, 2*i + 1, MPFR_RNDN);
    mpfr_div_si(x, x, 2, MPFR_RNDN);
    mpfr_mul(x, x, dx, MPFR_RNDN); /* x = (i + 1/2) * dx */
    mpfr_mul(arr[i], in[i], x, MPFR_RNDN);
    /* arr[i] = in[i] * x; */
  }

  /* do the sine transform */
  mpsint11(arr, npt, NULL);

  /* 1. the following dx is from integration
   * 2. divided by dk is for 1/k */
  /* fac = fac0 * dx / dk */
  mpfr_mul(fac, fac0, dx, MPFR_RNDN);
  mpfr_div(fac, fac, dk, MPFR_RNDN);

  for ( i = 0; i < npt; i++ ) { /* form out(k) / k */
    mpfr_set_si(x, 2*i + 1, MPFR_RNDN);
    mpfr_div_si(x, x, 2, MPFR_RNDN);
    mpfr_div(x, fac, x, MPFR_RNDN); /* fac/(i + 1/2) */
    mpfr_mul(out[i], arr[i], x, MPFR_RNDN);
    /* out[i] = arr[i] * fac / i */
  }
  mpfr_clear(fac);
  mpfr_clear(x);
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, const char *sdr, int doHNC)
{
  mpfr_t dr, dk, B2, facr2k, fack2r, tmp1, tmp2;
  mpfr_t Bc0, dBc, Bv0, dBv, eps;
  mpfr_t **ck, **tk, *fr, *fk, *tl, *cl, *Bc, *Bv;
  mpfr_t *powtlr, **yr0, **yr, *arr, *vc;
  int i, dm, l, u, j, jl, k;

  mpfr_init(dr);
  mpfr_init(dk);
  mpfr_init(B2);
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

  /* the sign transform already includes a factor of 2
   * so we only need to multiply 2*pi for the forward transform
   * in the backward case, we need to divide 2*pi by (2*pi)^dim */
  mpfr_const_pi(facr2k, MPFR_RNDN);
  mpfr_mul_si(facr2k, facr2k, 2, MPFR_RNDN); /* facr2k = 2*pi */

  mpfr_sqr(fack2r, facr2k, MPFR_RNDN);
  mpfr_si_div(fack2r, 1, fack2r, MPFR_RNDN); /* fack2r = 1/(2*pi)^2 */

  MAKE1DARR(arr, npt);

  MAKE1DARR(fr, npt);
  MAKE1DARR(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    mpfr_set_si(fr[i], (i < dm) ? -1 : 0, MPFR_RNDN);
  sphr11(npt, fr, fk, dr, dk, facr2k, arr); /* f(r) --> f(k) */

  /* B2 = 2*PI/3; */
  mpfr_const_pi(B2, MPFR_RNDN);
  mpfr_mul_si(B2, B2, 2, MPFR_RNDN);
  mpfr_div_si(B2, B2, 3, MPFR_RNDN);

  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  mpfr_set(Bv[2], B2, MPFR_RNDN);
  mpfr_set(Bc[2], B2, MPFR_RNDN);

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ )
    mpfr_set(ck[0][i], fk[i], MPFR_RNDN);

  MAKE1DARR(tl, npt);
  MAKE1DARR(cl, npt);

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
      mpfr_set_si(tl[i], 0, MPFR_RNDN); /* tl[i] = 0; */
    for ( u = 0; u < l; u++ )
      for ( i = 0; i < npt; i++ ) {
        /* tl[i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]); */
        mpfr_add(tmp1, ck[u][i], tk[u][i], MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, ck[l-1-u][i], MPFR_RNDN);
        mpfr_add(tl[i], tl[i], tmp1, MPFR_RNDN);
      }
    for ( i = 0; i < npt; i++ )
      mpfr_set(tk[l][i], tl[i], MPFR_RNDN); /* tk[l][i] = tl[i] */

    /* t_l(k) --> t_l(r) */
    sphr11(npt, tl, tl, dk, dr, fack2r, arr);

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
        /* powtlr = tl^j/j! */
        for ( i = 0; i < npt; i++) {
          mpfr_div_si(tmp1, tl[i], l, MPFR_RNDN);
          mpfr_mul(powtlr[i], powtlr[i], tmp1, MPFR_RNDN);
          /* powtlr[i] *= tl[i]/l */
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
        mpfr_neg(cl[i], tl[i], MPFR_RNDN); /* cl[i] = -tl[i]; */
      for ( i = dm; i < npt; i++ ) {
        mpfr_sub(cl[i], yr[l][i], tl[i], MPFR_RNDN);
        /* cl[i] = yr[l][i] - tl[i]; */
      }

      if ( mkcorr && l >= 2 ) { /* try to correct */
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
          mpfr_set(vc[i], cl[i], MPFR_RNDN);
        }

        mpfr_set_si(Bc0, 0, MPFR_RNDN);
        mpfr_set_si(dBc, 0, MPFR_RNDN);
        for ( i = 0; i < npt; i++ ) {
          mpfr_set_si(tmp1, 2*i + 1, MPFR_RNDN);
          mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
          mpfr_mul(tmp1, tmp1, dr, MPFR_RNDN); /* tmp1 = (i + 1/2) * dr  */

          mpfr_sqr(tmp1, tmp1, MPFR_RNDN);
          mpfr_mul(tmp2, cl[i], tmp1, MPFR_RNDN);
          mpfr_add(Bc0, Bc0, tmp2, MPFR_RNDN);
          /* Bc0 += cl[l] * ((i + 1/2) * dr)^2; */

          mpfr_add_si(tmp2, fr[i], 1, MPFR_RNDN);
          mpfr_mul(tmp2, tmp2, tmp1, MPFR_RNDN);
          mpfr_mul(tmp2, tmp2, vc[i], MPFR_RNDN);
          mpfr_add(dBc, dBc, tmp2, MPFR_RNDN);
          /* dBc += vc[i] * (1 + fr[i]) * ((i + 1/2) * dr)^2 */
        }
        mpfr_const_pi(tmp1, MPFR_RNDN);
        mpfr_mul_si(tmp1, tmp1, 4, MPFR_RNDN);
        mpfr_mul(tmp1, tmp1, dr, MPFR_RNDN);
        mpfr_div_si(tmp1, tmp1, l + 2, MPFR_RNDN);
        mpfr_neg(tmp1, tmp1, MPFR_RNDN);
        mpfr_mul(Bc0, Bc0, tmp1, MPFR_RNDN);
        mpfr_mul(dBc, dBc, tmp1, MPFR_RNDN);
        /*
        Bc0 *= -4*PI*dr/(l+2);
        dBc *= -4*PI*dr/(l+2);
        */
        mpfr_mul_si(tmp1, yr[l][dm], 2, MPFR_RNDN);
        mpfr_sub(tmp1, tmp1, yr[l][dm+1], MPFR_RNDN);
        mpfr_mul(Bv0, B2, tmp1, MPFR_RNDN);
        mpfr_mul_si(tmp1, vc[dm], 2, MPFR_RNDN);
        mpfr_sub(tmp1, tmp1, vc[dm+1], MPFR_RNDN);
        mpfr_mul(dBv, B2, tmp1, MPFR_RNDN);
        /*
        Bv0 = B2*(2*yr[l][dm] - yr[l][dm+1]);
        dBv = B2*(2*vc[dm] - vc[dm+1]);
        */

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
          mpfr_add(cl[i], cl[i], tmp2, MPFR_RNDN);
          /* cl[i] += (fr[i] + 1) * eps * vc[i] */
        }
      }
      /* Bv[l+2] = B2*(yr[l][dm] + yr[l][dm+1])/2; */
      mpfr_add(tmp1, yr[l][dm], yr[l][dm - 1], MPFR_RNDN);
      mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
      /* Bv[l+2] = B2*(2*yr[l][dm] - yr[l][dm+1]); */
      /*
      mpfr_mul_si(tmp1, yr[l][dm], 2, MPFR_RNDN);
      mpfr_sub(tmp1, tmp1, yr[l][dm + 1], MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
      */

    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ )
        mpfr_neg(cl[i], tl[i], MPFR_RNDN); /* cl[i] = -tl[i]; */
      for ( i = dm; i < npt; i++ )
        mpfr_set_si(cl[i], 0, MPFR_RNDN); /* cl[i] = 0; */
      /* Bv(l+2) = B2 * tl(1+) */
      /* Bv[l+2] = B2*(tl[dm] + tl[dm-1])/2; */
      mpfr_add(tmp1, tl[dm], tl[dm-1], MPFR_RNDN);
      mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
      /* Bv[l+2] = B2*(2*tl[dm] - tl[dm+1]); */
      /*
      mpfr_mul_si(tmp1, tl[dm], 2, MPFR_RNDN);
      mpfr_sub(tmp1, tmp1, tl[dm + 1], MPFR_RNDN);
      mpfr_mul(Bv[l+2], B2, tmp1, MPFR_RNDN);
      */
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    mpfr_set_si(Bc[l+2], 0, MPFR_RNDN);
    for ( i = 0; i < npt; i++ ) {
      mpfr_set_si(tmp1, i*2 + 1, MPFR_RNDN);
      mpfr_div_si(tmp1, tmp1, 2, MPFR_RNDN);
      mpfr_mul(tmp1, tmp1, dr, MPFR_RNDN); /* (i + 1/2) * dr */
      mpfr_sqr(tmp1, tmp1, MPFR_RNDN);
      mpfr_mul(tmp1, tmp1, cl[i], MPFR_RNDN);
      mpfr_add(Bc[l+2], Bc[l+2], tmp1, MPFR_RNDN);
      /* Bc[l+2] += cl[i] * ((i + 1/2) * dr)^2; */
    }
    /* Bc[l+2] *= -4*PI*dr/(l+2); */
    mpfr_const_pi(tmp1, MPFR_RNDN);
    mpfr_mul_si(tmp1, tmp1, -4, MPFR_RNDN);
    mpfr_mul(tmp1, tmp1, dr, MPFR_RNDN);
    mpfr_div_si(tmp1, tmp1, l + 2, MPFR_RNDN);
    mpfr_mul(Bc[l+2], Bc[l+2], tmp1, MPFR_RNDN);

    /* tmp1 = Bc/B2^(l+1), tmp2 = Bv/B2^(l+1) */
    mpfr_pow_si(tmp1, B2, l+1, MPFR_RNDN);
    mpfr_div(tmp2, Bv[l+2], tmp1, MPFR_RNDN);
    mpfr_div(tmp1, Bc[l+2], tmp1, MPFR_RNDN);
    mpfr_printf("Bc(%3d) = %20.10Re (%20.12Re), "
                "Bv(%3d) = %20.10Re (%20.12Re)\n",
                l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c(r) --> c(k) */
    sphr11(npt, cl, cl, dr, dk, facr2k, arr);

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ )
      mpfr_set(ck[l][i], cl[i], MPFR_RNDN); /* ck[l][i] = cl[i]; */
  }
  FREE1DARR(arr, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(fk, npt);
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
  FREE1DARR(tl, npt);
  FREE1DARR(cl, npt);
  FREE1DARR(Bc, nmax + 1);
  FREE1DARR(Bv, nmax + 1);

  mpfr_clear(dr);
  mpfr_clear(dk);
  mpfr_clear(B2);
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
