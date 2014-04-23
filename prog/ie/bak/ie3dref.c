/* Computing the virial coefficients of the 3D hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ie3dmp.c -lmpfr
 * */
#include <stdio.h>
#include <math.h>
#include "fft.h"
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"


#define PI 3.1415926535897932384626433832795
int nmax = 10;
char *dr = NULL;
int numpt = 1024;
int doHNC = 0;
int mkcorr = 0;
int type = 1;



#define MAKE1DARR(arr, n) { int i_; \
  xnew(arr, n); \
  for (i_ = 0; i_ < (n); i_++) \
    arr[i_] = 0; }

#define FREE1DARR(arr, n) { int i_; \
  for (i_ = 0; i_ < (n); i_++) \
    arr[i_] = 0; \
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
  argopt_add(ao, "-b", NULL, &dr, "interval of r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if (dr == NULL) dr = "0.01";
  argopt_dump(ao);
  argopt_close(ao);
}



/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr00(int npt, double *in, double *out,
    double dx, double dk, double fac0, double *arr)
{
  int i;
  double fac, x;

  for ( i = 1; i < npt; i++ ) { /* form in(x) * x */
    x = dx * i;
    arr[i] = in[i] * x;
  }

  /* do the sine transform */
  sint00(arr, npt);

  /* 1. the following dx is from integration
   * 2. divided by dk is for 1/k */
  fac = fac0 * dx / dk;

  for ( i = 1; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / i;
}



/* compute
 *    out(k) = 2*fac0/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr11(int npt, double *in, double *out,
    double dx, double dk, double fac0, double *arr)
{
  int i;
  double fac, x;

  for ( i = 0; i < npt; i++ ) { /* form in(x) * x */
    x = dx * (i + .5);
    arr[i] = in[i] * x;
  }

  /* do the sine transform */
  sint11(arr, npt);

  /* 1. the following dx is from integration
   * 2. divided by dk is for 1/k */
  fac = fac0 * dx / dk;

  for ( i = 0; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / (i + .5);
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, const char *sdr, int doHNC)
{
  double dr, dk, B2, facr2k, fack2r, tmp1, tmp2;
  double Bc0, dBc, Bv0, dBv, eps, idel;
  double **ck, **tk, *fr, *fk, *tl, *cl, *Bc, *Bv;
  double **yr, **yr0, **tr, *powtl, *arr, *vc;
  int i, im, dm, l, u, j, jl, k;

  dr = atof(sdr);
  dk = PI/npt/dr;
  im = (type == 0) ? 1 : 0;
  idel = (type == 0) ? 0 : 0.5;

  /* the bin index at which the hard core stops */
  dm = (int) (1/atof(sdr) + .5);
  printf("dr %g, dk %g, dm %d\n", dr, dk, dm);

  /* the sign transform already includes a factor of 2
   * so we only need to multiply 2*pi for the forward transform
   * in the backward case, we need to divide 2*pi by (2*pi)^dim */
  facr2k = 2 * PI;
  fack2r = 1 / (4 * PI * PI);

  MAKE1DARR(arr, npt);

  MAKE1DARR(fr, npt);
  MAKE1DARR(fk, npt);
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    fr[i] = (i < dm) ? -1 : 0;
  if (type == 0) {
    sphr00(npt, fr, fk, dr, dk, facr2k, arr); /* f(r) --> f(k) */
  } else {
    sphr11(npt, fr, fk, dr, dk, facr2k, arr); /* f(r) --> f(k) */
  }

  B2 = 2*PI/3;

  MAKE1DARR(Bc, nmax + 1);
  MAKE1DARR(Bv, nmax + 1);
  Bv[2] = B2;
  Bc[2] = B2;

  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  for ( i = 0; i < npt; i++ )
    ck[0][i] = fk[i];

  MAKE1DARR(tl, npt);
  MAKE1DARR(cl, npt);

  if ( doHNC ) {
    MAKE1DARR(powtl, npt);
    MAKE2DARR(yr, nmax + 1, npt);
    MAKE2DARR(yr0, nmax + 1, npt);
    for ( i = 0; i < npt; i++ ) yr0[0][i] = yr[0][i] = 1;
    if ( mkcorr ) {
      MAKE2DARR(tr, nmax + 1, npt);
      MAKE1DARR(vc, npt);
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
    if (type == 0) {
      sphr00(npt, tl, tl, dk, dr, fack2r, arr);
    } else {
      sphr11(npt, tl, tl, dk, dr, fack2r, arr);
    }

    if ( doHNC ) {
      if ( mkcorr )
        for ( i = 0; i < npt; i++ )
          tr[l][i] = tl[i];

      /* Hypernetted chain approximation
       * y(r) = exp(t(r))
       * c(r) = (1 + f(r)) y(r) - t(r) - 1 */
      /* back up, only need to save yr[1..nmax-l] */
      for ( u = 0; u < nmax - l; u++ ) {
        for ( i = 0; i < npt; i++ ) {
          yr0[u][i] = yr[u][i];
        }
      }
      //memcpy(yr0[1], yr[1], sizeof(mydouble) * (nmax - l) * npt);
      /* yr = yr0 * exp(rho^l t_l)
       *    = (yr0_0 + yr0_1 rho + yr0_2 rho^2 + ... )
       *    * (1 + t_l rho^l + t_{2l} rho^(2l)/2! + ... )
       * y[l...n](r) are updated */
      for ( i = 0; i < npt; i++ ) {
        powtl[i] = 1;
        /* powtl[i] = 1; */
      }
      for ( j = 1; j*l <= nmax; j++ ) {
        /* powtl = tl^j/j! */
        for ( i = 0; i < npt; i++) {
          powtl[i] *= tl[i]/l;
        }
        for ( jl = j * l, k = 0; k + jl <= nmax; k++ )
          /* yr_{k + jl} +=  yr0_k * tl^j/j! */
          for ( i = 0; i < npt; i++ ) {
            yr[jl+k][i] += yr0[k][i] * powtl[i];
          }
      }

      /* compute c(r) = (1 + f) y(r) - (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) {
        cl[i] = -tl[i];
      }
      for ( i = dm; i < npt; i++ ) {
        cl[i] = yr[l][i] - tl[i];
      }

      if ( mkcorr && l >= 2 ) { /* try to correct */
        for ( i = 0; i < npt; i++ ) {
          vc[i] = cl[i];
          /* vc[i] = tr[l][i] - yr[l][i]; */
        }

        Bc0 = 0;
        dBc = 0;
        for ( i = im; i < npt; i++ ) {
          tmp1 = dr * (i + idel);
          tmp1 = tmp1 * tmp1;
          Bc0 += cl[l] * tmp1;

          tmp2 = fr[i] + 1;
          tmp2 *= tmp1 * vc[i];
          dBc += tmp2;
          /* dBc += vc[i] * (1 + fr[i]) * (i * dr)^2 */
        }
        Bc0 *= -4*PI*dr/(l+2);
        dBc *= -4*PI*dr/(l+2);
        if (type == 0) {
          Bv0 = B2*(2*yr[l][dm] - yr[l][dm+1]);
          dBv = B2*(2*vc[dm] - vc[dm+1]);
        } else {
          Bv0 = B2*(yr[l][dm] + yr[l][dm-1])/2;
          dBv = B2*(vc[dm] + vc[dm-1])/2;
        }

        eps = -(Bv0 - Bc0) / (dBv - dBc);
        for ( i = im; i < npt; i++ ) {
          yr[l][i] += eps * vc[i];
          cl[i] += eps * vc[i] * (fr[i] + 1);
        }
      }
      if (type == 0) {
        Bv[l+2] = B2*(2*yr[l][dm] - yr[l][dm+1]);
      } else {
        Bv[l+2] = B2*(yr[l][dm] + yr[l][dm-1])/2;
      }

    } else {
      /* Percus-Yevick approximation
       * c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < dm; i++ ) {
        cl[i] = -tl[i];
      }
      for ( i = dm; i < npt; i++ ) {
        cl[i] = 0;
      }
      //Bv[l+2] = B2*tl[dm-1];
      if (type == 0) {
        Bv[l+2] = B2*(2*tl[dm] - tl[dm+1]);
      } else {
        Bv[l+2] = B2*(tl[dm] + tl[dm-1])/2;
      }
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    Bc[l+2] = 0;
    for ( i = im; i < npt; i++ ) {
      tmp1 = dr * (i + idel);
      tmp1 *= tmp1;
      Bc[l+2] += cl[i] * tmp1;
    }
    Bc[l+2] *= -4*PI*dr/(l+2);
    tmp1 = pow(B2, l+1);
    tmp2 = Bv[l+2]/tmp1;
    tmp1 = Bc[l+2]/tmp1;
    printf("Bc(%3d) = %20.10e (%20.12e), "
           "Bv(%3d) = %20.10e (%20.12e)\n",
            l+2, Bc[l+2], tmp1, l+2, Bv[l+2], tmp2);
    /* c(r) --> c(k) */
    if (type == 0) {
      sphr00(npt, cl, cl, dr, dk, facr2k, arr);
    } else {
      sphr11(npt, cl, cl, dr, dk, facr2k, arr);
    }

    /* save c(k) for the following iterations */
    for ( i = 0; i < npt; i++ ) {
      ck[l][i] = cl[i];
    }
  }
  FREE1DARR(arr, npt);
  FREE1DARR(fr, npt);
  FREE1DARR(fk, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(ck, nmax - 1, npt);
  if ( doHNC ) {
    FREE1DARR(powtl, npt);
    FREE2DARR(yr, nmax + 1, npt);
    FREE2DARR(yr0, nmax + 1, npt);
    if ( mkcorr ) {
      FREE1DARR(vc, npt);
      FREE2DARR(tr, nmax + 1, npt);
    }
  }
  FREE1DARR(tl, npt);
  FREE1DARR(cl, npt);
  FREE1DARR(Bc, nmax + 1);
  FREE1DARR(Bv, nmax + 1);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, dr, doHNC);
  return 0;
}
