/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * by the PY or HNC integral equations
 *  gcc iegsl.c -lgsl -lgslcblas
 * This program works for both even and odd dimensions
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"



#include "slowdht.h"
/*
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>
#include "xdouble.h"
*/

#include "ieutil.h"



#ifdef D
int dim = D;
#else
int dim = 2;
#endif

int nmax = 12;
double rmax = 0;
xdouble Rmax = 0;
int numpt = 1024;
int dohnc = 0;
int singer = 0;
int ring = 0;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;
int dhtdisk = SLOWDHT_USEDISK;
int snapshot = 0;

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
char systitle[32];

int ietype = 0;

xdouble hncamp = 0, hncq = 1; /* Marucho-Pettitt */
xdouble hncalpha = -1; /* Rogers-Young */
xdouble hcs = 1; /* Hutchinson-Conkie s */
xdouble rowphi = 0; /* Rowlinson's Phi */
xdouble invphi = 0; /* inverse Rowlinson's Phi */
xdouble verleta = 0, verletb = 0; /* Verlet modified */
xdouble bbpgs = 15./8; /* MS/BBPG s */

xdouble shift = 0, shiftinc = 0;
int shiftl0 = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  xdouble shiftn = 0;
  int usems = 0;

  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (exact)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-T", "%d", &ietype, "type of closure, 0: PY, 1: HNC-like; set the respective parameters automatically sets the type");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted-chain (HNC) like approximation");
  argopt_add(ao, "-q", "%" XDBLSCNF "f", &hncq,  "q of the hypernetted-chain (Marucho-Pettitt) approximation, y(r) = a0 exp(q t(r))");
  argopt_add(ao, "-a", "%" XDBLSCNF "f", &hncamp, "a0 of the hypernetted-chain (Marucho-Pettitt) approximation, 0: to be set by 1/q");
  argopt_add(ao, "--rya", "%" XDBLSCNF "f", &hncalpha, "alpha, in the Roger-Young switch function [1 - exp(-alpha*r)]^m");
  argopt_add(ao, "--hcs", "%" XDBLSCNF "f", &hcs, "s, in the Hutchinson-Conkie closure, y(r) = (1 + s t(r))^(1/s)");
  argopt_add(ao, "--rphi", "%" XDBLSCNF "f", &rowphi, "phi of the Rowlinson approximation, t(r) = (1 - phi) (y(r) - 1) + phi ln y(r)");
  argopt_add(ao, "--iphi", "%" XDBLSCNF "f", &invphi, "phi of the inverse Rowlinson approximation, y(r) = exp(t(r)) phi + (1 - phi) (1 + t(r))");
  argopt_add(ao, "--va", "%" XDBLSCNF "f", &verleta, "A of the Verlet approximation, y(r) = exp[t(r) - A t(r)^2 /2 / (1 + B t(r) / 2) ]");
  argopt_add(ao, "--vb", "%" XDBLSCNF "f", &verletb, "B of the Verlet approximation, y(r) = exp[t(r) - A t(r)^2 /2 / (1 + B t(r) / 2) ]");
  argopt_add(ao, "--bbpgs", "%" XDBLSCNF "f", &bbpgs, "s of the BBPG approximation, y(r) = exp[ (1 + s t(r))^(1/s) - 1 ]");
  argopt_add(ao, "--ms", "%b", &usems, "Martynov-Sarkisov approximation, y(r) = exp[ (1 + 2 t(r))^(1/2) - 1 ]");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "-c", "%" XDBLSCNF "f", &shift, "shift of t(r) in computing Bv");
  argopt_add(ao, "-C", "%" XDBLSCNF "f", &shiftn, "shift of t(r) in computing Bv, in terms of n");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &shiftinc, "increment of the shift");
  argopt_add(ao, "-L", "%d", &shiftl0, "minimal order l for the shift");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "--disk", "%d", &dhtdisk, "use disk for discrete Hankel transform (DHT)");
  argopt_add(ao, "--tmpdir", "%s", &slowdht_tmpdir, "directory for temporary files");
  argopt_add(ao, "--dhtblock", "%u", &slowdht_block, "block size for DHT input/output");
  argopt_add(ao, "--dhtquiet", "%b", &slowdht_quiet, "suppress DHT output");
  argopt_add(ao, "-s", "%b", &snapshot, "save intermediate snapshots");
  argopt_add(ao, "-G", "%b", &gaussf, "Gaussian model instead of hard spheres");
  argopt_add(ao, "-I", "%d", &invexp, "exponent of the inverse potential r^(-n)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);

  if ( rmax <= 0 ) rmax = nmax + 2;

  /* closure type */
  if ( hncq > 0 && hncamp <= 0 ) /* Marucho-Pettitt extended HNC */
    hncamp = 1/hncq; /* set amp automatically = 1/q */
  if ( argopt_isset(ao, hncamp) || argopt_isset(ao, hncq)
    || argopt_isset(ao, hncalpha) ) {
    ietype = IETYPE_HNC;
    dohnc = 1; /* turn on HNC, if necessary */
  }

  if ( argopt_isset(ao, hcs) ) ietype = IETYPE_HC;
  if ( argopt_isset(ao, rowphi) ) ietype = IETYPE_ROWLINSON;
  if ( argopt_isset(ao, invphi) ) ietype = IETYPE_INVROWLINSON;
  if ( argopt_isset(ao, verleta) || argopt_isset(ao, verletb) ) ietype = IETYPE_VERLET;
  if ( argopt_isset(ao, usems) ) { ietype = IETYPE_BBPG; bbpgs = 2; }
  if ( argopt_isset(ao, bbpgs) ) ietype = IETYPE_BBPG;

  if ( ietype > IETYPE_HNC ) dohnc = 0;

  if ( dohnc ) ietype = IETYPE_HNC;

  if ( mkcorr || ietype > IETYPE_HNC ) /* Singer and ring formulas are inapplicable to corrections */
    singer = ring = 0;
  if ( singer ) ring = 1;

  /* correction parameters */
  if ( mkcorr && ietype == IETYPE_PY ) { dohnc = 1; }
  if ( argopt_isset(ao, shiftn) ) shift = shiftn + shiftinc*2;
  /* for the Gaussian model, Bv is exact for l <= 2 (n <= 4)
   * otherwise, Bv is exact for l = 1 (n = 3) */
  if ( shiftl0 <= 0 ) shiftl0 = gaussf ? 3 : 2;

  /* models */
  if ( gaussf || invexp > 0 ) smoothpot = 1;
  if ( gaussf ) strcpy(systitle, "hGF");
  else if ( invexp ) sprintf(systitle, "hINV%d", invexp);
  else strcpy(systitle, "h");

  printf("D %d, rmax %f, npt %d, ietype %d, HNC %d, a0 %g, q %g, %s\n",
      dim, (double) rmax, numpt, ietype, dohnc,
      (double) hncamp, (double) hncq, systitle);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute
 *  out(k) = fac Int {from 0 to infinity} dr
 *           r^(D/2) in(r) J_{D/2-1}(k r) */
static void sphr(xdouble *in, xdouble *out, xdouble fac,
    xdht *dht, xdouble *arr, xdouble *r2p, xdouble *k2p)
{
  int i, npt = dht->size;

  for ( i = 0; i < npt; i++ ) arr[i] = in[i] * r2p[i];
  XDHT(apply)(dht, arr, out);
  for ( i = 0; i < npt; i++ ) out[i] *= fac / k2p[i];
}



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, int dohnc)
{
  xdouble facr2k, fack2r, surfr, surfk;
  xdouble Bc, Bv, Bm = 0, Bh = 0, Br = 0, B2, tmp1, tmp2, fcorr = 0;
  xdouble *fr, *rdfr = NULL, *swr = NULL;
  xdouble *crl, *trl, **cr = NULL, **tr = NULL, **ck, **tk;
  xdouble **yr = NULL, *yrl = NULL, *arr, *vc = NULL, *yrcoef = NULL;
  xdouble *ri, *ki, *r2p, *k2p, *rDm1, *kDm1;
  int i, dm, l, l0 = 1;
  xdht *dht;
  clock_t t0 = clock(), t1;

  dht = XDHT(newx)(npt, (xdouble) dim/2 - 1, rmax, dhtdisk);

  for ( tmp1 = 0, dm = 0; dm < npt; dm++, tmp1 = tmp2 )
    if ((tmp2 = XDHT(x_sample)(dht, dm)) > 1) {
      printf("r %g, %g, dm %d. \n", (double) tmp1, (double) tmp2, dm);
      break;
    }

  facr2k = POW(PI*2, (xdouble) dim/2);
  /* the factor (kmax/xmax)^2 is adapted for the inverse
   * discrete hankel transform using gsl_dht */
  fack2r = pow_si(dht->kmax / dht->xmax, 2);
  fack2r /= POW(PI*2, (xdouble) dim/2);

  /* B2 = (1/2) (PI*2)^(dim/2) / dim!! for an even dim
   *    = (PI*2)^((dim - 1)/2) / dim!! for an odd dim */
  B2 = dim % 2 ? 1 : 0.5;
  for ( i = 2 + dim % 2; i <= dim; i += 2 ) B2 *= PI*2/i;
  surfr = B2 * 2 * dim;
  tmp2 = dht->kmax / dht->xmax;
  surfk = surfr * tmp2 * tmp2 / pow_si(PI*2, dim);
  if ( gaussf ) B2 = SQRT(pow_si(PI, dim))/2;
  /* TODO: compute B2 for inverse potential */
  if ( invexp > 0 ) B2 = 0;
  printf("D %d, dm %d, rmax %f, ietype %d, B2 %.10e\n",
      dim, dm, (double) rmax, ietype, (double) B2);

  MAKE1DARR(ri,   npt);
  MAKE1DARR(ki,   npt);
  MAKE1DARR(r2p,  npt);
  MAKE1DARR(k2p,  npt);
  MAKE1DARR(rDm1, npt);
  MAKE1DARR(kDm1, npt);
  for ( i = 0; i < npt; i++ ) {
    ri[i]   = XDHT(x_sample)(dht, i);
    r2p[i]  = POW(ri[i], (xdouble) dim/2 - 1);
    ki[i]   = XDHT(k_sample)(dht, i);
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

  /* auxiliary array for the Hankel transform */
  MAKE1DARR(arr,  npt);

  MAKE1DARR(fr, npt);
  if ( smoothpot ) MAKE1DARR(rdfr, npt);
  MAKE1DARR(swr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE1DARR(yrl, npt);
  MAKE2DARR(ck, nmax - 1, npt);
  MAKE2DARR(tk, nmax - 1, npt);

  /* compute f(r) and f(k) = c0(k) */
  for ( i = 0; i < npt; i++ ) { /* compute f(r) = exp(-beta u(r)) - 1 */
    if ( gaussf ) { /* f(r) = exp(-r^2) */
      xdouble r2 = ri[i] * ri[i];
      fr[i] = -EXP(-r2);
      rdfr[i] = -2*r2*fr[i];
    } else if ( invexp > 0 ) { /* inverse potential r^(-invexp) */
      xdouble pot = POW(ri[i], -invexp);
      fr[i] = EXP(-pot) - 1;
      rdfr[i] = pot * invexp * (fr[i] + 1);
    } else { /* hard-sphere */
      fr[i] = (i < dm) ? -1 : 0;
    }
    crl[i] = fr[i];
  }

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt);
  }

  if ( ietype > IETYPE_HNC || singer ) {
    MAKE2DARR(tr, nmax - 1, npt);
  }

  if ( ietype > IETYPE_HNC ) { /* coefficients */
    MAKE1DARR(yrcoef, nmax - 1);
    if ( ietype == IETYPE_HC ) {
      init_hccoef(yrcoef, nmax - 1, hcs);
    } else if ( ietype == IETYPE_BBPG ) {
      init_bbpgcoef(yrcoef, nmax - 1, bbpgs);
    } else if ( ietype == IETYPE_ROWLINSON ) {
      init_rowlinsoncoef(yrcoef, nmax - 1, rowphi);
    } else if ( ietype == IETYPE_INVROWLINSON ) {
      init_invrowlinsoncoef(yrcoef, nmax - 1, invphi);
    } else if ( ietype == IETYPE_VERLET ) {
      init_verletcoef(yrcoef, nmax - 1, verleta, verletb);
    }
    print_yrcoef(yrcoef, 7);
  }

  if ( dohnc ) { /* initialize the Rogers-Young switch function */
    for ( i = 0; i < npt; i++ ) {
      swr[i] = 1;
      if ( hncalpha >= 0 )
        swr[i] = 1 - EXP( -hncalpha * ri[i] );
    }
  }

  if ( dohnc || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) yr[0][i] = hncamp / swr[i];
  }

  if ( mkcorr ) {
    MAKE1DARR(vc, npt);
  }

  t1 = clock();

  if ( snapshot )
    l0 = snapshot_open(dim, nmax, rmax, dohnc, mkcorr, ring, singer,
        npt, ck, tk, cr, tr, crl, trl, yr);

  fnvir = savevirheadx(fnvir, systitle, dim, l0, nmax,
      ietype, mkcorr, npt, rmax, t1 - t0,
      hncamp, hncq, hncalpha, shift, shiftinc, shiftl0);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* c_l(r) --> c_l(k) for the previous l */
    sphr(crl, ck[l-1], facr2k, dht, arr, r2p, k2p);

    if ( ring ) {
      /* compute the ring sum based on ck */
      Bh = get_ksum(l, npt, ck, kDm1, &Br);
      Br = (dohnc ? -Br * (l+1) : -Br * 2) / l;
    }

    /* compute t_l(k) from c_0(k), ... c_{l-1}(k) */
    get_tk_oz(l, npt, ck, tk);

    /* t_l(k) --> t_l(r) */
    sphr(tk[l], trl, fack2r, dht, arr, k2p, r2p);

    if ( tr != NULL ) { /* save tr if needed */
      COPY1DARR(tr[l], trl, npt);
    }

    /* compute the cavity function y(r) */
    if ( dohnc ) { /* specialized HNC closure */
      get_yr_hncx(l, nmax, npt, yr, trl, hncq, swr);
      for ( i = 0; i < npt; i++ )
        yrl[i] = yr[l][i] + (1 - hncamp * hncq) * trl[i];
    } else if ( ietype == IETYPE_PY) { /* PY closure */
      for ( i = 0; i < npt; i++ )
        yrl[i] = trl[i];
    } else {
      get_yr_series(l, npt, tr, yrcoef, yrl);
    }

    for ( i = 0; i < npt; i++ ) /* c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      crl[i] = (fr[i] + 1) * (mkcorr ? trl[i] : yrl[i]) - trl[i];

    Bv = get_Bv(npt, mkcorr ? trl : yrl, smoothpot, rdfr, rDm1, dim, dm, B2);
    Bc = get_Bc(l, npt, crl, rDm1);

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( dohnc ) {
        Bm = get_Bm_singer(l, npt, cr, tr, rDm1);
        Bh = get_Bh_singer(l, npt, cr, tr, rDm1) - Bh*(l+1)/2;
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, rDm1);
        Bh = get_Bp_py(l, npt, cr, tr, rDm1) - Bh;
      }
    } else {
      Bm = Bh = 0;
    }

    if ( mkcorr ) {
      /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yrl[i] - trl[i];

      /* apply the shift */
      if ( l >= shiftl0 ) Bv *= 1 + shift + shiftinc * l;

      Bm = get_corr1x(l, npt, dm, crl, fr, rdfr, rDm1,
                      dim, B2, vc, &Bc, &Bv, &fcorr);
    }

    savevir(fnvir, dim, l+2, Bc, Bv, Bm, Bh, Br, B2, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, ri, crl, trl, vc, yr);
    if ( snapshot )
      snapshot_take(l, npt, ck[l-1], tk[l], crl, trl, nmax, yr);
  }
  savevirtail(fnvir, clock() - t1);

  FREE1DARR(arr,  npt);
  FREE1DARR(ri,   npt);
  FREE1DARR(ki,   npt);
  FREE1DARR(r2p,  npt);
  FREE1DARR(k2p,  npt);
  FREE1DARR(rDm1, npt);
  FREE1DARR(kDm1, npt);
  FREE1DARR(fr,   npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(swr,  npt);
  FREE1DARR(crl,  npt);
  FREE1DARR(trl,  npt);
  FREE1DARR(yrl,  npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(vc, npt);
  FREE1DARR(yrcoef, nmax - 1);
  XDHT(free)(dht);
  return 0;
}



/* adjust rmax such that r = 1 lies at the middle of the dm'th and dm+1'th bins */
static xdouble jadjustrmax(double rmax0, int npt)
{
  xdouble dr, km, kp, kM, rmax;
  double nu = dim*.5 - 1;
  int dm;

  dr = rmax0 / npt;
  dm = (int)(1/dr + .5);
  kM = gsl_sf_bessel_zero_Jnu(nu, npt + 1);
  while ( 1 ) {
    km = gsl_sf_bessel_zero_Jnu(nu, dm);
    kp = gsl_sf_bessel_zero_Jnu(nu, dm + 1);
    /* adjust rmax such that rmax (j_{nu,k} * .5 + j_{nu,k+1} * .5) / j_{nu,M} = 1 */
    rmax = kM*2/(km + kp);
    //printf("dm %d, rmax %g, k %g, %g, %g\n", dm, (double) rmax, (double) km, (double) kp, (double) kM);
    if ( rmax >= rmax0 - 1e-8 || dm == 1 ) break;
    dm--;
  }
  dr = rmax / npt;
  printf("D %d, %d bins, %d within the hard core (dr %g), rmax %g, k %g - %g\n",
      dim, npt, dm, (double) dr, (double) rmax, (double) km, (double) kp);
  return rmax;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  if ( Rmax <= 0 ) Rmax = jadjustrmax(rmax, numpt);
  intgeq(nmax, numpt, Rmax, dohnc);
  return 0;
}
