/* Computing the virial coefficients of the 5D hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ie5d.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ie5d.c -lfftw3l
 * */
#include <stdio.h>
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
  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 5D hard-sphere fluid";
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-b", "%" DBLSCNF "f", &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  printf("rmax = %" DBLPRNF "f, HNC %d\n", dr * numpt, doHNC);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*(2*pi)^2 Int {from 0 to infinity} in(r) r^4
 *          [sin(k r)/(k r)^3 - cos(k r)/(k r)^2] dr
 * if forward = 1.
 * otherwise compute
 *    out(r) = 2/(2*pi)^3 Int {from 0 to infinity} in(k) k^4
 *          [sin(k r)/(k r)^3 - cos(k r)/(k r)^2] dk
 * `arr' are used in the intermediate steps by the FFTW plans `p'
 * */
static void sphr(int npt, mydouble *in, mydouble *out, mydouble dx,
    int forward, FFTWPFX(plan) p[2], mydouble *arr)
{
  int i;
  mydouble fac, r, k, dk;

  /* for the sine transform */
  for ( i = 0; i < npt; i++ ) { /* * x */
    r = dx*(2*i + 1)/2;
    arr[i] = in[i] * r;
  }
  FFTWPFX(execute)(p[0]);
  /* 1. the sign transform already includes a factor of 2
   *    so we only need to multiply 2*pi for the forward transform
   *    in the backward case, we need to divide 2*pi by (2*pi)^dim
   * 2. the following dx is from integration
   * 3. the final pi/npt/dx is dk which is to be used for 1/k */
  dk = PI/npt/dx;
  fac = (forward ? pow(2*PI, 2) : 1/pow(2*PI, 3)) * dx;
  for ( i = 0; i < npt; i++ ) { /* / k^3 */
    k = dk*(2*i + 1)/2;
    out[i] = arr[i] * fac / (k * k * k);
  }

  /* for the cosine transform */
  for ( i = 0; i < npt; i++ ) { /* x^2 */
    r = dx*(2*i + 1)/2;
    arr[i] = -in[i] * r * r;
  }
  FFTWPFX(execute)(p[1]);
  for ( i = 0; i < npt; i++ ) { /* / k^2 */
    k = dk*(2*i + 1)/2;
    out[i] += arr[i] * fac / (k * k);
  }
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, mydouble dr, int doHNC)
{
  mydouble dk = PI/dr/npt, B2, tmp1, tmp2;
  mydouble **ck, **tk, *fk, *fr, *Bc, *Bv;
  mydouble *arr, *tlk, *tlr, *clk, *clr;
  mydouble *powtlr, **yr0, **yr;
  FFTWPFX(plan) plans[2];
  int i, dm, l, u, j, jl, k;

  dm = (int) (1/dr + .5);

  /* fftw auxiliary array */
  xnew(arr, npt);
  /* sine transform */
  plans[0] = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);
  /* cosine transform */
  plans[1] = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_REDFT11, FFTW_ESTIMATE);

  xnew(fr, npt);
  xnew(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    fr[i] = (i < dm) ? -1. : 0;
  sphr(npt, fr, fk, dr, 1, plans, arr); /* f(r) --> f(k) */

  B2 = pow(2*PI, 2)/(3*5);
  xnew(Bc, nmax + 1);
  xnew(Bv, nmax + 1);
  Bc[2] = Bv[2] = B2;

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ ) {
    ck[0][i] = fk[i];
    tk[0][i] = 0;
  }

  xnew(clk, npt);
  xnew(clr, npt);
  xnew(tlk, npt);
  xnew(tlr, npt);

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
    sphr(npt, tlk, tlr, dk, 0, plans, arr);

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
      Bv[l+2] = pow(2*PI, 2)/(3*5)*(1.5*yr[l][dm] - .5*yr[l][dm+1]);
    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) clr[i] = -tlr[i];
      for ( i = dm; i < npt; i++ ) clr[i] = 0;
      Bv[l+2] = pow(2*PI, 2)/(3*5)*(1.5*tlr[dm] - .5*tlr[dm+1]);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 2 (2 pi)^K/(2*K - 1)!! r^(2*K) dr */
    for ( Bc[l+2] = 0, i = 0; i < npt; i++ ) {
      mydouble r = dr * (2*i + 1)/2;
      Bc[l+2] += clr[i] * (r * r * r * r);
    }
    Bc[l+2] *= -2*pow(2*PI, 2)/3 * dr/(l+2);
    tmp1 = pow(B2, -l-1);
    tmp2 = Bv[l+2] * tmp1;
    tmp1 = Bc[l+2] * tmp1;
    printf("Bc(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e), "
           "Bv(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e)\n",
           l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c_l(r) --> c_l(k) */
    sphr(npt, clr, clk, dr, 1, plans, arr);

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ ) ck[l][i] = clk[i];
  }
  free(arr);
  FFTWPFX(destroy_plan)(plans[0]);
  FFTWPFX(destroy_plan)(plans[1]);
  free(clr); free(clk);
  free(tlr); free(tlk);
  free(fr); free(fk);
  FREE2DARR(tk); FREE2DARR(ck);
  if ( doHNC ) {
    free(powtlr); FREE2DARR(yr); FREE2DARR(yr0);
  }
  free(Bc); free(Bv);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, dr, doHNC);
  FFTWPFX(cleanup)();
  return 0;
}
