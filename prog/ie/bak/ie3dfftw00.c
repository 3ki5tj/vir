/* Computing the virial coefficients of the 3D hard-sphere fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ie3d.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ie3d.c -lfftw3l
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

int mkcorr = 0;



#define MAKE2DARR(arr, n1, n2) { int l_; \
  xnew(arr, n1); xnew(arr[0], (n1) * (n2)); \
  for ( l_ = 1; l_ < (n1); l_++ ) arr[l_] = arr[0] + l_ * (n2); }

#define FREE2DARR(arr) { free(arr[0]); free(arr); }



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-b", "%" DBLSCNF "f", &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  printf("rmax = %" DBLPRNF "f, HNC %d\n", dr * numpt, doHNC);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, mydouble *in, mydouble *out, mydouble dx,
    mydouble fac0, FFTWPFX(plan) p, mydouble *arr)
{
  int i;
  mydouble fac;

  for ( i = 1; i < npt; i++ ) /* form in(x) * x */
    arr[i] = in[i] * dx * i;
  FFTWPFX(execute)(p);
  /* 1. the following dx is from integration
   * 2. the final pi/npt/dx is dk which is to be used for 1/k */
  fac = fac0 * dx / (PI/npt/dx);
  for ( i = 1; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / i;
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, mydouble dr, int doHNC)
{
  mydouble dk = PI/dr/npt, B2, facr2k, fack2r, tmp1, tmp2;
  mydouble **ck, **tk, *fr, *fk, *tl, *cl, *Bc, *Bv;
  mydouble **yr, **yr0, **tr, *powtl, *arr, *vc;
  int i, dm, l, u, j, jl, k;
  FFTWPFX(plan) plan;

  dm = (int) (1/dr + .5);

  /* the sign transform already includes a factor of 2
   * so we only need to multiply 2*pi for the forward transform
   * in the backward case, we need to divide 2*pi by (2*pi)^dim */
  facr2k = 2*PI;
  fack2r = 1/(4*PI*PI);
  xnew(arr, npt);
  plan = FFTWPFX(plan_r2r_1d)(npt - 1, arr + 1, arr + 1, FFTW_RODFT00, FFTW_ESTIMATE);

  xnew(fr, npt);
  xnew(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    fr[i] = (i < dm) ? -1. : 0;
  sphr(npt, fr, fk, dr, facr2k, plan, arr); /* f(r) --> f(k) */

  B2 = 2*PI/3;
  xnew(Bc, nmax + 1);
  xnew(Bv, nmax + 1);
  Bc[2] = Bv[2] = B2;

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ ) {
    ck[0][i] = fk[i];
    tk[0][i] = 0;
  }

  xnew(tl, npt);
  xnew(cl, npt);

  if ( doHNC ) {
    xnew(powtl, npt);
    MAKE2DARR(yr, nmax + 1, npt);
    MAKE2DARR(yr0, nmax + 1, npt);
    for ( i = 0; i < npt; i++ ) yr0[0][i] = yr[0][i] = 1;
    for ( j = 1; j <= nmax; j++ ) /* coefficient of rho^j */
      for ( i = 0; i < npt; i++ ) yr[j][i] = 0;
    if ( mkcorr ) {
      MAKE2DARR(tr, nmax + 1, npt);
      for ( j = 0; j <= nmax; j++ )
        for ( i = 0; i < npt; i++ ) tr[j][i] = 0;
      xnew(vc, npt);
    }
  }

  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute t_l(k) */
    for ( i = 0; i < npt; i++ ) tl[i] = 0;
    for ( u = 0; u < l; u++ )
      for ( i = 0; i < npt; i++ )
        tl[i] += ck[l-1-u][i] * (ck[u][i] + tk[u][i]);
    for ( i = 0; i < npt; i++ )
      tk[l][i] = tl[i];

    /* t_l(k) --> t_l(r) */
    sphr(npt, tl, tl, dk, fack2r, plan, arr);

    if ( doHNC ) {
      if ( mkcorr )
        for ( i = 0; i < npt; i++ ) tr[l][i] = tl[i];

      /* Hypernetted chain approximation
       * y(r) = exp(t(r))
       * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
      /* back up, only need to save yr[1..nmax-l] */
      memcpy(yr0[1], yr[1], sizeof(mydouble) * (nmax - l) * npt);
      /* yr = yr0 * exp(rho^l t_l)
       *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
       *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
       * y[l...n](r) are updated */
      for ( i = 0; i < npt; i++ ) powtl[i] = 1;
      for ( j = 1; j*l <= nmax; j++ ) {
        /* powtl = tl^j/j! */
        for ( i = 0; i < npt; i++) powtl[i] *= tl[i]/j;
        for ( jl = j * l, k = 0; k + jl <= nmax; k++ )
          /* yr_{k + jl} +=  yr0_k * tl^j/j! */
          for ( i = 0; i < npt; i++ )
            yr[jl+k][i] += yr0[k][i] * powtl[i];
      }

      /* compute c(r) = (1 + f) y(r) - (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) cl[i] = -tl[i];
      for ( i = dm; i < npt; i++ ) cl[i] = yr[l][i] - tl[i];

      if ( mkcorr && l >= 2 ) { /* try to correct */
        mydouble Bc0 = 0, dBc = 0, Bv0 = 0, dBv = 0, eps;

        for ( i = 0; i < npt; i++ )
          vc[i] = cl[i];
          /* vc[i] = tr[l][i] - yr[l][i]; */

        for ( i = 0; i < npt; i++ ) {
          Bc0 += cl[i] * i * i;
          dBc += vc[i] * (1 + fr[i]) * i * i;
        }
        Bc0 *= -4*PI*(dr*dr*dr)/(l+2);
        dBc *= -4*PI*(dr*dr*dr)/(l+2);
        Bv0 = B2*(2*yr[l][dm] - yr[l][dm+1]);
        dBv = B2*(2*vc[dm] - vc[dm+1]);
        eps = -(Bv0 - Bc0) / (dBv - dBc);
        for ( i = 0; i < npt; i++ ) {
          cl[i] += eps * vc[i] * (1 + fr[i]);
          yr[l][i] += eps * vc[i];
        }
      }
      Bv[l+2] = B2*(2*yr[l][dm] - yr[l][dm+1]);
    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) cl[i] = -tl[i];
      for ( i = dm; i < npt; i++ ) cl[i] = 0;
      Bv[l+2] = B2*(2*tl[dm] - tl[dm+1]);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr
     * In the 3D case the current formula appeared to work better than:
       Bc[l+2] += cl[i] * (i*(i+1)+1./3); */
    for ( Bc[l+2] = 0, i = 0; i < npt; i++ )
      Bc[l+2] += cl[i] * i * i;
    Bc[l+2] *= -4*PI*(dr*dr*dr)/(l+2);
    tmp1 = pow(B2, -l-1);
    tmp2 = Bv[l+2] * tmp1;
    tmp1 = Bc[l+2] * tmp1;
    printf("Bc(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e), "
           "Bv(%3d) = %20.10" DBLPRNF "e (%20.12" DBLPRNF "e)\n",
           l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c(r) --> c(k) */
    sphr(npt, cl, cl, dr, facr2k, plan, arr);

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ ) ck[l][i] = cl[i];
  }
  FFTWPFX(destroy_plan)(plan); free(arr);
  free(fr); free(fk); FREE2DARR(tk); FREE2DARR(ck);
  if ( doHNC ) {
    free(powtl); FREE2DARR(yr); FREE2DARR(yr0);
    if ( mkcorr ) {
      free(vc); FREE2DARR(tr);
    }
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
