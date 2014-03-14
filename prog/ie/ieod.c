/* Computing the virial coefficients of a odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ieod.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ieod.c -lfftw3l
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#ifdef LDBL
#define PI 3.1415926535897932384626433832795L
typedef long double mydouble;
#define FFTWPFX(f) fftwl_##f
#define DBLSCNF "L"
#define DBLPRNF "L"
#else
#define PI 3.1415926535897932384626433832795
typedef double mydouble;
#define FFTWPFX(f) fftw_##f
#define DBLSCNF "l"
#define DBLPRNF ""
#endif



int D = 3;
int K = 1; /* (D - 1)/2 */

#ifdef LDBL
int nmax = 54; /* maximal order */
mydouble dr = (mydouble) 0.0002L; /* grid point */
int numpt = 1<<19; /* # of points along r */
#else
int nmax = 10;
mydouble dr = 0.01;
int numpt = 1024;
#endif
int doHNC = 0;



#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); xnew(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr) { free(arr[0]); free(arr); }



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &D, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-b", "%" DBLSCNF "f", &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if (D < 3 || D % 2 == 0) argopt_help(ao);
  K = (D - 1)/2;
  printf("D %d, K %d, rmax %" DBLPRNF "f, HNC %d\n", D, K, dr * numpt, doHNC);
  argopt_close(ao);
}




/* get the coefficient c_l of spherical Bessel function jn(x)
 *  jn(x) = Sum_{l = 0 to n} c_l [sin(x) or cos(x)] / x^{n + l + 1},
 * where the l'th term is sin(x) if n + l is even, or cos(x) otherwise */
static void getjn(int *c, int n)
{
  int i, k;
  char *fs[2] = {"sin(x)", "cos(x)"};

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


/* compute
 *    out(k) = 2 fac0 Int {from 0 to infinity} dr
 *             r^(2K) in(r) j_{D-1}(k r)/(k r)^{D - 1}
 * `arr' are used in the intermediate steps by the FFTW plans `p'
 * */
static void sphr(int npt, mydouble *in, mydouble *out, mydouble dx,
    mydouble fac0, FFTWPFX(plan) p[2], mydouble *arr[2], int *coef,
    mydouble **r2p, mydouble **k2q)
{
  int i, l, iscos;
  mydouble fac = fac0*dx, dk = PI/npt/dx;

  /* clear the output */
  for ( i = 0; i < npt; i++ ) out[i] = 0;

  /* several rounds of transforms */
  for ( l = 0; l < K; l++ ) {
    /* decide if we want to do a sine or cosine transform */
    iscos = (K + l + 1) % 2;
    for ( i = 0; i < npt; i++ ) {
      //arr[iscos][i] = in[i] * pow(dx*(2*i + 1)/2, K - l) * coef[l];
      arr[iscos][i] = in[i] * coef[l] * r2p[l][i];
    }
    FFTWPFX(execute)(p[iscos]);
    for ( i = 0; i < npt; i++ ) {
      //out[i] += arr[iscos][i] * fac / pow(dk*(2*i + 1)/2, K + l);
      out[i] += arr[iscos][i] * fac * k2q[l][i];
    }
  }
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, mydouble dr, int doHNC)
{
  mydouble dk = PI/dr/npt, B2;
  mydouble **ck, **tk, *fk, *fr, *Bc, *Bv;
  int i, dm, l, u, j, jl, k;
  mydouble *arrs[2], *tlk, *tlr, *clk, *clr;
  FFTWPFX(plan) plans[2];
  mydouble **yr, **yr0, *powtlr;
  mydouble **r2p, **invr2p, **k2p, **invk2p, *r2Dm1;
  mydouble facr2k, fack2r;
  int *coef;
  dm = (int) (1/dr + .5);

  /* compute the coefficients of the spherical Bessel function */
  xnew(coef, K);
  getjn(coef, K - 1);

  /* FFTW auxiliary array */
  xnew(arrs[0], npt); /* sine transform */
  plans[0] = FFTWPFX(plan_r2r_1d)(npt, arrs[0], arrs[0], FFTW_RODFT11, FFTW_ESTIMATE);
  xnew(arrs[1], npt); /* cosine transform */
  plans[1] = FFTWPFX(plan_r2r_1d)(npt, arrs[1], arrs[1], FFTW_REDFT11, FFTW_ESTIMATE);

  MAKE2DARR(r2p, K, npt) /* r^{K - l} for l = 0, ..., K */
  MAKE2DARR(k2p, K, npt) /* k^{K - l} for l = 0, ..., K */
  MAKE2DARR(invr2p, K, npt) /* r^{-K-l} for l = 0, ..., K */
  MAKE2DARR(invk2p, K, npt) /* k^{-K-l} for l = 0, ..., K */
  for ( i = 0; i < npt; i++ ) {
    mydouble r = dr * (2*i + 1)/2, rl = 1;
    mydouble k = dk * (2*i + 1)/2, kl = 1;
    mydouble invr = 1/r, invrl, invk = 1/k, invkl;
    for ( l = 1; l <= K; l++ ) {
      r2p[K - l][i] = (rl *= r);
      k2p[K - l][i] = (kl *= k);
    }
    invrl = 1/rl;
    invkl = 1/kl;
    for ( l = 0; l < K; l++, invrl *= invr, invkl *= invk ) {
      invr2p[l][i] = invrl;
      invk2p[l][i] = invkl;
    }
  }

  /* compute r^(D - 1) dr */
  xnew(r2Dm1, npt);
  for (i = 0; i < npt; i++) {
    mydouble r = dr * (2*i + 1)/2;
    r2Dm1[i] = pow(r, D - 1) * dr;
  }

  /* B2 = (2*PI)^K/(2 K + 1)!! */
  B2 = 1;
  facr2k = 1;
  for (i = 1; i <= K; i++) {
    B2 *= 2*PI/(2*i + 1);
    facr2k *= 2*PI;
  }
  fack2r = 1/(facr2k * 2*PI);
  xnew(Bc, nmax + 1);
  xnew(Bv, nmax + 1);
  Bc[2] = Bv[2] = B2;

  /* construct f(r) and f(k) */
  xnew(fr, npt);
  xnew(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    fr[i] = (i < dm) ? -1. : 0;
  sphr(npt, fr, fk, dr, facr2k, plans, arrs, coef, r2p, invk2p); /* f(r) --> f(k) */

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ ) {
    ck[0][i] = fk[i];
    tk[0][i] = 0;
  }

  xnew(clk, npt); xnew(clr, npt);
  xnew(tlk, npt); xnew(tlr, npt);

  if ( doHNC ) {
    xnew(powtlr, npt);
    MAKE2DARR(yr, nmax + 1, npt);
    MAKE2DARR(yr0, nmax + 1, npt);
    for ( i = 0; i < npt; i++ ) yr0[0][i] = yr[0][i] = 1;
    for ( j = 1; j <= nmax; j++ ) /* coefficient of rho^j */
      for ( i = 0; i < npt; i++ ) yr[j][i] = 0;
  }

  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute t_l(k) */
    for ( i = 0; i < npt; i++ ) tlk[i] = 0;
    for ( u = 0; u < l; u++ )
      for ( i = 0; i < npt; i++ )
        tlk[i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]);
    for ( i = 0; i < npt; i++ )
      tk[l][i] = tlk[i];

    /* t_l(k) --> t_l(r) */
    sphr(npt, tlk, tlr, dk, fack2r, plans, arrs, coef, k2p, invr2p);

    if ( doHNC ) {
      /* Hypernetted chain approximation
       * y(r) = exp(t(r))
       * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
      /* back up, only need to save yr[1..nmax-l] */
      memcpy(yr0[1], yr[1], sizeof(mydouble) * (nmax - l) * npt);
      /* yr = yr0 * exp(rho^l t_l)
       *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
       *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
       * y[l...n](r) are updated */
      for ( i = 0; i < npt; i++ ) powtlr[i] = 1;
      for ( j = 1; j*l <= nmax; j++ ) {
        /* powtl(r) = tl(r)^j/j! */
        for ( i = 0; i < npt; i++) powtlr[i] *= tlr[i]/j;
        for ( jl = j * l, k = 0; k + jl <= nmax; k++ )
          /* yr_{k + jl} +=  yr0_k * tl^j/j! */
          for ( i = 0; i < npt; i++ )
            yr[jl+k][i] += yr0[k][i] * powtlr[i];
      }

      /* compute c(r) = (1 + f) y(r) - (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) clr[i] = -tlr[i];
      for ( i = dm; i < npt; i++ ) clr[i] = yr[l][i] - tlr[i];
      Bv[l+2] = B2*(3*yr[l][dm] - yr[l][dm+1])/2;
    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) clr[i] = -tlr[i];
      for ( i = dm; i < npt; i++ ) clr[i] = 0;
      Bv[l+2] = B2*(3*tlr[dm] - tlr[dm+1])/2;
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 2 (2 pi)^K/(2*K - 1)!! r^(2*K) dr
     *           = -2*B2*D/(l+2) Int c_l(r) r^(2*K) dr */
    for ( Bc[l+2] = 0, i = 0; i < npt; i++ ) {
      Bc[l+2] += clr[i] * r2Dm1[i];
    }
    Bc[l+2] *= -2*B2*D/(l+2);
    printf("Bc(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e), "
           "Bv(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e)\n",
           l+2, Bc[l+2], Bc[l+2]/pow(B2, l+1),
           l+2, Bv[l+2], Bv[l+2]/pow(B2, l+1));
    /* c_l(r) --> c_l(k) */
    sphr(npt, clr, clk, dr, facr2k, plans, arrs, coef, r2p, invk2p);

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ ) ck[l][i] = clk[i];
  }
  FFTWPFX(destroy_plan)(plans[0]); free(arrs[0]);
  FFTWPFX(destroy_plan)(plans[1]); free(arrs[1]);
  free(clr); free(clk);
  free(tlr); free(tlk);
  free(fr); free(fk);
  FREE2DARR(tk); FREE2DARR(ck);
  if ( doHNC ) {
    free(powtlr); FREE2DARR(yr); FREE2DARR(yr0);
  }
  free(Bc); free(Bv);
  free(coef);
  free(r2Dm1);
  FREE2DARR(r2p); FREE2DARR(invr2p);
  FREE2DARR(k2p); FREE2DARR(invk2p);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, dr, doHNC);
  FFTWPFX(cleanup)();
  return 0;
}
