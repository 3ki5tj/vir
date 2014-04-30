/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc ieedgsl.c -lgsl -lgslcblas
 * This program works for both even and odd dimensions
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"


typedef double xdouble;
#include "ieutil.h"



#ifdef D
int dim = D;
#else
int dim = 2;
#endif

int nmax = 10;
double rmax = 10.24;
int numpt = 512;
int doHNC = 0;
int singer = 0;
int mkcorr = 0;
int verbose = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &dim, "dimension (odd) integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", "%lf", &rmax, "rmax");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &doHNC, "use the hypernetted chain approximation");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}




/* compute
 *  out(k) = fac Int {from 0 to infinity} dr
 *           r^(D/2) in(r) J_{D/2-1}(k r) */
static void sphr(double *in, double *out, double fac,
    gsl_dht *dht, double *arr, double *r2p, double *k2p)
{
  int i, npt = dht->size;

  for ( i = 0; i < npt; i++ ) arr[i] = in[i] * r2p[i];
  gsl_dht_apply(dht, arr, out);
  for ( i = 0; i < npt; i++ ) out[i] *= fac / k2p[i];
}



/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, double rmax, int doHNC)
{
  double facr2k, fack2r, surfr, surfk;
  double Bc, Bv, Bm, Bh, Br, B2, B2p, tmp1, tmp2;
  double *fr, *crl, *trl, **ck, **tk, **cr = NULL, **tr = NULL;
  double **yr = NULL, *arr, *vc = NULL;
  double *r2p, *k2p, *rDm1, *kDm1;
  int i, dm, l;
  gsl_dht *dht;

  dht = gsl_dht_new(npt, dim*.5 - 1, rmax);

  for ( tmp1 = 0, dm = 0; dm < (int) dht->size; dm++, tmp1 = tmp2 )
    if ((tmp2 = gsl_dht_x_sample(dht, dm)) > 1) {
      printf("r %g, %g, dm %d\n", tmp1, tmp2, dm);
      break;
    }

  /* auxiliary array */
  MAKE1DARR(arr, npt);

  MAKE1DARR(r2p, dht->size + 1)
  MAKE1DARR(k2p, dht->size + 1)
  for ( i = 0; i < (int) dht->size; i++ ) {
    r2p[i] = pow(gsl_dht_x_sample(dht, i), dim*.5 - 1);
    k2p[i] = pow(gsl_dht_k_sample(dht, i), dim*.5 - 1);
  }

  facr2k = pow(PI*2, dim*.5);
  /* the factor (kmax/xmax)^2 is adapted for the inverse
   * discrete hankel transform using gsl_dht */
  fack2r = pow(dht->kmax / dht->xmax, 2);
  fack2r /= pow(PI*2, dim*.5);

  /* B2 = (1/2) (PI*2)^(dim/2) / dim!! for an even dim
   *    = (PI*2)^((dim - 1)/2) / dim!! for an odd dim */
  B2 = dim % 2 ? 1 : 0.5;
  for ( i = 2 + dim % 2; i <= dim; i+=2 ) B2 *= PI*2/i;
  surfr = B2 * 2 * dim;
  tmp2 = dht->kmax / dht->xmax;
  surfk = surfr * tmp2 * tmp2 / pow(PI*2, dim);
  printf("B2 %g, r2k %g, k2r %g\n", B2, facr2k, fack2r);

  /* compute r^(dim - 1) dr, used for integration
   *   \int r dr / xmax^2
   * ==>
   *   2/(j_{nu,M})^2 * 1/[J_{nu+1}(j_{nu,k})]^2
   * r = j_{nu,k} / j_{nu,M} */
  MAKE1DARR(rDm1, npt);
  MAKE1DARR(kDm1, npt);
  for ( i = 0; i < npt; i++ ) {
    tmp1 = gsl_dht_x_sample(dht, i);
    tmp2 = pow(tmp1, dim - 2) * surfr;
    rDm1[i] = tmp2 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
    tmp1 = gsl_dht_k_sample(dht, i);
    tmp2 = pow(tmp1, dim - 2) * surfk;
    kDm1[i] = tmp2 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
  }

  MAKE1DARR(fr, npt);
  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);

  /* construct f(r) and f(k) */
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    fr[i] = (i < dm) ? -1. : 0;
  sphr(fr, ck[0], facr2k, dht, arr, r2p, k2p); /* f(r) --> f(k) */

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt);
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( doHNC || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) yr[0][i] = 1;
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  B2p = B2;
  for ( l = 1; l < nmax - 1; l++ ) {
    /* compute the ring sum based on ck */
    Bh = get_ksum(l, npt, ck, kDm1, &Br);
    Br = (doHNC ? -Br * (l+1) : -Br * 2) / l;

    /* compute t_l(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(tk[l], trl, fack2r, dht, arr, k2p, r2p);

    if ( tr != NULL ) COPY1DARR(tr[l], trl, npt);

    if ( yr != NULL ) { /* compute the cavity function y(r) */
      get_yr_hnc(l, nmax, npt, yr, trl);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yr[l][i] - trl[i];
    }

    if ( doHNC ) {
      /* HNC approximation: c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = (fr[i] + 1) * yr[l][i] - trl[i];
      if ( mkcorr ) {
        for ( i = 0; i < npt; i++ )
          vc[i] = yr[l][i] - trl[i];
      }
      Bv = contactv(yr[l], dm, B2);
    } else {
      /* PY approximation: c(r) = f(r) (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = fr[i] * trl[i];
      Bv = contactv(trl, dm, B2);
    }

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( doHNC ) {
        Bm = get_Bm_singer(l, npt, cr, tr, rDm1);
        Bh = get_Bh_singer(l, npt, cr, tr, rDm1) - Bh*(l+1)/2;
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, rDm1);
        Bh = get_Bp_py(l, npt, cr, tr, rDm1) - Bh;
      }
    }

    if ( mkcorr ) {
      get_corr1_hnc_hs(l, npt, dm, doHNC ? yr[l] : trl, crl, fr, rDm1, B2, vc);
      Bv += contactv(vc, dm, B2);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) S_D r^(D-1) dr */
    Bc = -integr(npt, crl, rDm1) / (l + 2);
    B2p *= B2;
    printf("Bc(%3d) = %16.9e (%16.9e), "
           "Bv(%3d) = %16.9e (%16.9e), "
           "Bm(%3d) = %16.9e (%16.9e)\n",
           l+2, Bc, Bc/B2p, l+2, Bv, Bv/B2p, l+2, Bm, Bm/B2p);
    printf("Bh(%3d) = %16.9e (%16.9e), "
           "Br(%3d) = %16.9e (%16.9e)\n",
           l+2, Bh, Bh/B2p, l+2, Br, Br/B2p);

    /* c_l(r) --> c_l(k) */
    sphr(crl, ck[l], facr2k, dht, arr, r2p, k2p);
  }

  FREE1DARR(arr, npt);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE1DARR(fr, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(r2p, npt);
  FREE1DARR(k2p, npt);
  if ( yr != NULL ) FREE2DARR(yr, nmax - 1, npt);
  if ( vc != NULL ) FREE1DARR(vc, npt);
  gsl_dht_free(dht);
  return 0;
}



/* adjust rmax such that r = 1 lies at the middle of the dm'th and dm+1'th bins */
static void adjustrmax(double *rmax, int npt)
{
  double dr, km, kp, kM, nu = dim*.5 - 1;
  int dm;

  dr = *rmax / npt;
  dm = (int)(1/dr + .5);
  km = gsl_sf_bessel_zero_Jnu(nu, dm);
  kp = gsl_sf_bessel_zero_Jnu(nu, dm + 1);
  kM = gsl_sf_bessel_zero_Jnu(nu, npt + 1);
  /* adjust rmax such that rmax (j_{nu,k} * .5 + j_{nu,k+1} * .5) / j_{nu,M} = 1 */
  *rmax = kM*2/(km + kp);
  printf("%d bins (dr %g) within the hard core, rmax %g, k %g - %g\n",
      dm, dr, *rmax, km, kp);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  adjustrmax(&rmax, numpt);
  intgeq(nmax, numpt, rmax, doHNC);
  return 0;
}
