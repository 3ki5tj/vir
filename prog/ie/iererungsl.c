/* re-run iegsl.c from snapshots
 * This tools help fixing missing lines in hBnXXX.dat */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#ifndef DHT
#define DHT 1
#endif
#include "slowdht.h"
#include "ieutil.h"



#ifdef D
int dim = D;
#else
int dim = 2;
#endif

int nmax = 10;
double rmax = 0;
int numpt = 1024;
int dohnc = 0;
int singer = 1; /* always turn on singer */
int ring = 1; /* always turn on ring */
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "fixing GSL";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-R", "%lf", &rmax, "rmax");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted chain approximation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax <= 0 ) rmax = nmax + 2;
  if ( mkcorr ) singer = ring = 0;
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}




/* correct the hypernetted chain approximation
 * `vc' is the trial correction function to y(r)
 * The corresponding correction to c(r) is (f(r) + 1) vc(r)
 * which is only effective for r > 1 */
__inline static xdouble get_invcorr1_hs(int l, int npt, int dm,
    xdouble *yr, xdouble *cr, xdouble *fr, xdouble *rDm1,
    xdouble B2, xdouble *vc, xdouble *Bc0, xdouble *Bv0, xdouble *eps)
{
  int i;
  xdouble B, dBc, dBv;

  /* B = - Int c(r) dr / (l + 2) */
  B = dBc = 0;
  for ( i = 0; i < npt; i++ ) {
    B += cr[i] * rDm1[i];
    dBc += vc[i] * (1 + fr[i]) * rDm1[i];
  }
  B /= -(l + 2);
  dBc /= -(l + 2);
  *Bv0 = contactv(yr, dm, B2);
  dBv = contactv(vc, dm, B2);
  if ( l <= 1 ) {
    *Bc0 = B;
    *eps = 0;
    return *Bv0;
  }

  *eps = (B - (*Bv0)) / dBv;
  for ( i = 0; i < npt; i++ ) {
    vc[i] *= *eps;
    cr[i] -= (fr[i] + 1) * vc[i];
  }
  /* Sometimes, the *Bc0 and *eps values are wrong,
   * while B and *Bv0 are still correct.
   * This means cr[] and tr[] (in the PY case) are correct
   * but vc[], hence dBv and/or dBc, is incorrect */
  *Bc0 = ((dBv - dBc) * B + dBc * (*Bv0)) / dBv;
  return B;
}



static int rerun(int nmax, int npt, xdouble rmax, int dohnc)
{
  xdouble facr2k, fack2r, surfr, surfk;
  xdouble Bc = 0, Bv = 0, Bm = 0, Bh = 0, Br = 0, By = 0, B2, tmp1, tmp2, fcorr = 0;
  xdouble *fr, **cr, **tr, **ck, **tk, *yrl = NULL, *vc;
  xdouble *ri, *ki, *r2p, *k2p, *rDm1, *kDm1;
  xdouble *yr0, *lnyr0;
  int i, dm, l, l0 = 1;
  slowdht *dht;
  clock_t t0 = clock(), t1;

  dht = slowdht_newx(npt, (xdouble) dim/2 - 1, rmax,
      SLOWDHT_NOJJJ); /* we do not need Jjj */

  for ( tmp1 = 0, dm = 0; dm < npt; dm++, tmp1 = tmp2 )
    if ((tmp2 = slowdht_x_sample(dht, dm)) > 1) {
      printf("r %g, %g, dm %d. \n", (double) tmp1, (double) tmp2, dm);
      break;
    }

  facr2k = POW(PI*2, (xdouble) dim/2);
  fack2r = pow_si(dht->kmax / dht->xmax, 2);
  fack2r /= POW(PI*2, (xdouble) dim/2);

  B2 = dim % 2 ? 1 : 0.5;
  for ( i = 2 + dim % 2; i <= dim; i+=2 ) B2 *= PI*2/i;
  surfr = B2 * 2 * dim;
  tmp2 = dht->kmax / dht->xmax;
  surfk = surfr * tmp2 * tmp2 / pow_si(PI*2, dim);
  printf("D %d, B2 %g, r2k %g, k2r %g\n", dim, (double) B2, (double) facr2k, (double) fack2r);

  MAKE1DARR(ri,   npt);
  MAKE1DARR(ki,   npt);
  MAKE1DARR(r2p,  npt);
  MAKE1DARR(k2p,  npt);
  MAKE1DARR(rDm1, npt);
  MAKE1DARR(kDm1, npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i]   = slowdht_x_sample(dht, i);
    r2p[i]  = POW(ri[i], (xdouble) dim/2 - 1);
    ki[i]   = slowdht_k_sample(dht, i);
    k2p[i]  = POW(ki[i], (xdouble) dim/2 - 1);

    /* compute r^(dim - 1) dr, used for integration
     *   \int r dr / xmax^2
     * ==>
     *   2/(j_{nu,M})^2 * 1/[J_{nu+1}(j_{nu,k})]^2
     * r = j_{nu,k} / j_{nu,M} */
    tmp1 = pow_si(ri[i], dim - 2) * surfr;
    rDm1[i] = tmp1 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
    tmp1 = pow_si(ki[i], dim - 2) * surfk;
    kDm1[i] = tmp1 * 2 / (dht->kmax * dht->kmax * dht->J2[i+1]);
  }

  MAKE1DARR(yr0, nmax - 1);
  MAKE1DARR(lnyr0, nmax - 1);
  yr0[0] = 1;
  lnyr0[0] = 0;

  MAKE1DARR(fr, npt);
  MAKE2DARR(ck, nmax - 1, npt)
  MAKE2DARR(tk, nmax - 1, npt)
  MAKE2DARR(cr, nmax - 1, npt)
  MAKE2DARR(tr, nmax - 1, npt)
  MAKE1DARR(yrl, npt);
  MAKE1DARR(vc, npt);

  /* construct f(r) and f(k) */
  for ( i = 0; i < npt; i++ ) /* compute f(r) = exp(-beta u(r)) - 1 */
    cr[0][i] = fr[i] = (i < dm) ? -1. : 0;

  l0 = snapshot_open(dim, nmax, rmax, dohnc, mkcorr, ring, singer,
      npt, ck, tk, cr, tr, NULL, NULL, NULL);
  t1 = clock();
  fnvir = savevirhead(fnvir, "_h", dim, 1, nmax,
      dohnc, mkcorr, npt, rmax, t1 - t0);

  for ( l = 1; l < l0; l++ ) {
    /* compute the ring sum based on ck */
    if ( ring ) {
      Bh = get_BhBrk(l, npt, dohnc, ck, kDm1, ki, &Br);
    }

    /* compute the cavity function y(r) */
    get_yr_hnc(l, npt, yrl, tr);

    if ( mkcorr ) {
      /* vc(r) is unaffected */
      for ( i = 0; i < npt; i++ )
        vc[i] = yrl[i] - tr[l][i];

      /* Bv is unaffected */
      if ( dohnc ) {
        Bv = contactv(yrl, dm, B2);
        yr0[l] = get_zerosep(yrl, ri);
      } else {
        Bv = contactv(tr[l], dm, B2);
        yr0[l] = get_zerosep(tr[l], ri);
      }

      /* in the PY case, y(r) = 1 + t(r) */
      Bm = get_invcorr1_hs(l, npt, dm, dohnc ? yrl : tr[l],
          cr[l], fr, rDm1, B2, vc, &Bc, &Bv, &fcorr);
    } else {
      /* without correction */
      if ( dohnc ) {
        Bv = contactv(yrl, dm, B2);
        yr0[l] = get_zerosep(yrl, ri);
      } else {
        Bv = contactv(tr[l], dm, B2);
        yr0[l] = get_zerosep(tr[l], ri);
      }
      Bc = -integr(npt, cr[l], rDm1) / (l + 2);

      if ( dohnc ) {
        Bm = get_Bm_singer(l, npt, cr, tr, rDm1);
        Bh += get_Bh_singer(l, npt, cr, tr, rDm1);
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, rDm1);
      }
    }

    By = update_lnyr0(l, yr0, lnyr0);
    savevir(fnvir, dim, l+2, Bc, Bv, Bm, Bh, Br, By, B2, mkcorr, fcorr);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(ri,   npt);
  FREE1DARR(ki,   npt);
  FREE1DARR(r2p,  npt);
  FREE1DARR(k2p,  npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr,   npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE1DARR(yrl, npt);
  FREE1DARR(vc, npt);
  FREE1DARR(yr0, nmax - 1);
  FREE1DARR(lnyr0, nmax - 1);
  slowdht_free(dht);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  rerun(nmax, numpt, jadjustrmax(rmax, numpt, dim), dohnc);
  return 0;
}
