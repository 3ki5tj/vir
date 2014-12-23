/* Computing the virial coefficients of an odd-dimensional hard-sphere fluid
 * from the Yvon-Born-Green integral equations
 * Test compile (no external library)
 *  gcc ybgvir.c -lm
 * For the normal double precision
 *  gcc -DFFTW -DDHT ybgvir.c -lfftw3 -lgsl -lgslcblas
 * Or for the long double precision
 *  gcc -DFFTW -DDHT -DLDBL ybgvir.c -lfftw3l -lgsl -lgslcblas
 * Or for the 128-bit precision
 *  gcc -DFFTW -DDHT -DF128 ybgvir.c -lfftw3q -lgsl -lgslcblas -lquadmath -lm
 * */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



int dim = D;
int K; /* (dim - 1) / 2 for an odd dimension */
int nmax = 12;
xdouble T = 1;
xdouble beta = 1;
double rmax = 0;
xdouble Rmax = 0;
int numpt = 32768;
int ffttype = 1;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;
int snapshot = 0;

enum { CORR_HNC = 1, CORR_HC = 2,
  CORR_RHODEP = 3, CORR_RHODEPB4 = 4, CORR_RHODEPB5 = 5};

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
int dolj = 0; /* Lennard-Jones potential */
char systitle[32];

xdouble hcs = 0.56;


sphr_t *sphr;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the Yvon-Born-Green closure for hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (exact)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type, 0: integral grid points, 1: half-integer grid points");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "--corr", "%d", &mkcorr, "linearly correct the closure, 1: HNC, 2: HC, 3: Gaskell, 4: Gaskell B4, 5: Gaskell B4-B5");
  argopt_add(ao, "--hcs", "%" XDBLSCNF "f", &hcs, "s in the Hutchinson-Conkie closure");
#ifdef DHT
  argopt_add(ao, "--disk", "%d", &dhtdisk, "use disk for discrete Hankel transform (DHT)");
  argopt_add(ao, "--tmpdir", "%s", &slowdht_tmpdir, "directory for temporary files");
  argopt_add(ao, "--dhtblock", "%u", &slowdht_block, "block size for DHT input/output");
  argopt_add(ao, "--dhtquiet", "%b", &slowdht_quiet, "suppress DHT output");
#endif
  argopt_add(ao, "-s", "%b", &snapshot, "save intermediate snapshots");
  argopt_add(ao, "-G", "%b", &gaussf, "Gaussian model instead of hard spheres");
  argopt_add(ao, "-I", "%d", &invexp, "exponent of the inverse potential r^(-n)");
  argopt_add(ao, "-v", "%+", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);

  if ( dim < 1 ) argopt_help(ao);
  K = (dim - 1)/2;
  if ( dim % 2 == 0 ) {
    if ( !argopt_isset(ao, numpt) ) /* adjust the default number of points */
      numpt = 1024;
  }

  if ( rmax <= 0 ) rmax = nmax + 2;

  /* models */
  if ( gaussf || invexp > 0 ) smoothpot = 1;

  { /* write the name of the system */
    char syst[16] = "";
    if ( gaussf ) strcpy(syst, "GF");
    else if ( invexp ) sprintf(syst, "INV%d", invexp);
    sprintf(systitle, "%s%s", (dim % 2 ? "" : "h"), syst);
  }

  printf("D %d, Rmax %f, npt %d, %s\n",
      dim, (double) Rmax, numpt, systitle);

  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* compute w(k) from f(k), y(r = 1), and h(k) */
__inline static void get_wk_ybg(int l, int npt, xdouble *wkl,
    xdouble *fk, xdouble *yd, xdouble **hk)
{
  int i, u;

  for ( i = 0; i < npt; i++ )
    for ( wkl[i] = 0, u = 0; u < l; u++ )
      wkl[i] += fk[i] * yd[u] * hk[l - u - 1][i];
}



/* compute w(k) from f(k), y(r = 1), and h(k), self-consistent
 * with rho-dependent lambda */
__inline static void get_wk_ybgx(int l, int npt,
    xdouble *wkl, xdouble *lambda,
    xdouble *fk, xdouble *yd,
    xdouble **wk, xdouble **hk)
{
  int i, u, v;

  for ( i = 0; i < npt; i++ ) {
    for ( wkl[i] = 0, u = 0; u < l; u++ )
      wkl[i] += fk[i] * yd[u] * hk[l - u - 1][i];

    /* add low-order corrections from the closure */
    for ( v = 0; v < l - 2; v++ ) {
      xdouble y;
      for ( y = 0, u = 1; u < l - v; u++ )
        y += (hk[u][i] - wk[u][i] - fk[i]*yd[u]) * hk[l - v - u - 1][i];
      wkl[i] += y * lambda[v];
    }
  }
}



/* Gaskell 1968 Physics Letters, Vol 27A, number 4, pages 209-210 */
__inline static void mkcorrfunc0(int l, int npt,
    xdouble *vcr, xdouble *vck,
    xdouble *fk, xdouble *yd, xdouble **wk, xdouble **hk)
{
  int i;

  if (l < 2) return;
  for ( i = 0; i < npt; i++ )
    vck[i] = (hk[1][i] - wk[1][i] - fk[i]*yd[1]) * hk[0][i];
  sphr_k2r(sphr, vck, vcr); /* compute vcr */
}



__inline static void mkcorrfunc1(int l, int npt,
    xdouble *vcr, xdouble *vck, xdouble **wk, xdouble **hk)
{
  int i, u;

  for ( i = 0; i < npt; i++ )
    for ( vck[i] = -wk[l][i], u = 0; u < l; u++ )
      vck[i] += (hk[u][i] - wk[u][i]) * hk[l - u - 1][i];
  sphr_k2r(sphr, vck, vcr); /* compute vcr */
}



__inline static void mkcorrfunc2(int l, int npt,
    xdouble *vcr, xdouble *vck,
    xdouble **tr, xdouble *wrl, xdouble **hk, xdouble s)
{
  int i, u, v;
  xdouble *a, *b;

  xnew(a, l + 1);
  xnew(b, l + 1);
  for ( i = 0; i < npt; i++ ) {
    /* compute the indirect correlation function t(k)
     * from the Ornstein-Zernike relation */
    a[0] = 0;
    for ( v = 1; v <= l; v++ )
      for ( a[v] = 0, u = 0; u < v; u++ )
        a[v] += (hk[u][i] - a[u]) * hk[v-1-u][i];
    vck[i] = a[l];
  }
  sphr_k2r(sphr, vck, tr[l]); /* compute tr */
  for ( i = 0; i < npt; i++ ) {
    /* form 1 + t */
    for ( a[0] = 1, u = 1; u <= l; u++ )
      a[u] = s*tr[u][i];
    /* compute log(1 + s t)/s */
    log_series(l+1, a, b);
    vcr[i] = b[l]/s - wrl[i];
  }
  sphr_r2k(sphr, vcr, vck); /* vc(r) */
  free(a);
  free(b);
}



/* correct the YBG equation
 * this function is in fact identical to get_corr1_hs
 * with cr --> hr */
#define get_ybgcorr1_hs get_corr1_hs
/* y(r) = exp w(r), for the highest order correction
 * d y_l(r) = d w_l(r) = vc(r)
 * So
 *    dBv = contactv(vcr, dm, B2);
 */
/* d(beta*P)/drho = (1 - rho ck0) = 1/(1 + rho hk0)
 * for the highest order component
 * d [ (l+2) B_{l+2} rho^{l+1} ] = - d(hk0l) rho^{l+1}
 * d(hk0) = Int (1 + f) d(yrl) dr
 * So
 *    dBc = -integre(npt, vcr, fr, rDm1)/(l+2);
 */



/* compute virial coefficients from integral equations */
static int intgeq(int nmax, int npt, xdouble rmax, xdouble Rmax, int ffttype)
{
  xdouble Bc, Bv, Bm = 0, By = 0, B2, fcorr = 0;
  xdouble *fr, *fk, *rdfr = NULL;
  xdouble **wr, *hrl, *yrl, *yd;
  xdouble *vcr = NULL, *vck = NULL, **tr = NULL, *lambda = NULL;
  xdouble **wk = NULL, *wkl = NULL, **hk = NULL, *hk0;
  int i, l, l0 = 1;
  clock_t t0 = clock(), t1;

  (void) Rmax;
  sphr = sphr_open(dim, npt, rmax, Rmax, ffttype);
  B2 = gaussf ? sphr->B2g : sphr->B2hs;
  /* TODO: compute B2 for inverse potential */
  if ( invexp > 0 ) B2 = 0;
  printf("D %d, dr %f, dm %d, rmax %f, ffttype %d, B2 %.10e\n",
      dim, (double) sphr->rmax/npt, sphr->dm, (double) sphr->rmax,
      ffttype, (double) B2);

  MAKE1DARR(fr, npt);
  MAKE1DARR(fk, npt);
  if ( smoothpot ) MAKE1DARR(rdfr, npt);
  MAKE1DARR(yrl, npt);
  MAKE1DARR(hrl, npt);
  MAKE1DARR(wkl, npt);
  MAKE2DARR(wr, nmax - 1, npt);
  MAKE2DARR(hk, nmax - 1, npt);
  MAKE1DARR(yd, nmax - 1);
  MAKE1DARR(hk0, nmax - 1);

  yd[0] = 1;
  hk0[0] = -2*B2;

  /* compute f(r) */
  mkfr(npt, beta, NULL, NULL, NULL, fr, rdfr, sphr->ri, sphr->dm, gaussf, invexp, dolj);
  COPY1DARR(hrl, fr, npt);
  /* compute fk */
  sphr_r2k(sphr, fr, fk);
  COPY1DARR(hk[0], fk, npt); /* hk[0] = fk */

  if ( mkcorr ) {
    MAKE1DARR(vcr, npt);
    MAKE1DARR(vck, npt);
    MAKE2DARR(wk, nmax - 1, npt);
    if ( mkcorr == CORR_HC ) {
      MAKE2DARR(tr, nmax - 1, npt);
    }
  }
  MAKE1DARR(lambda, nmax - 1);

  t1 = clock();

  fnvir = savevirhead(fnvir, systitle, dim, l0, nmax,
      IETYPE_YBG, 0, npt, rmax, t1 - t0);

  for ( l = l0; l < nmax - 1; l++ ) {
    /* get w_l(k) */
    if ( mkcorr >= CORR_RHODEP ) {
      get_wk_ybgx(l, npt, wkl, lambda, fk, yd, wk, hk);
    } else {
      get_wk_ybg(l, npt, wkl, fk, yd, hk);
    }

    if ( wk != NULL ) {
      COPY1DARR(wk[l], wkl, npt);
    }

    /* w_l(k) --> w_l(r) */
    sphr_k2r(sphr, wkl, wr[l]);

    /* y_l(r) = [exp w(r)]_l */
    get_yr_hnc(l, npt, yrl, wr);

    for ( i = 0; i < npt; i++ ) {
      /* h(r) = (f(r) + 1) y(r) - 1 */
      hrl[i] = (fr[i] + 1) * yrl[i];
    }

    /* h_l(r) --> h_l(k) */
    sphr_r2k(sphr, hrl, hk[l]);

    Bv = get_Bv(npt, yrl, smoothpot, rdfr, sphr->rDm1, dim, sphr->dm, B2);
    Bc = get_Bc_hr(l, npt, hrl, hk0, sphr->rDm1);

    if ( mkcorr ) {
      if ( mkcorr == 1 ) {
        mkcorrfunc1(l, npt, vcr, vck, wk, hk);
      } else if ( mkcorr == 2 ) { /* HC */
        mkcorrfunc2(l, npt, vcr, vck, tr, wr[l], hk, hcs);
      } else { /* Gaskell */
        mkcorrfunc0(l, npt, vcr, vck, fk, yd, wk, hk);
      }

      if (  mkcorr <= CORR_RHODEP ||
           (mkcorr == CORR_RHODEPB4 && l <= 2) ||
           (mkcorr == CORR_RHODEPB5 && l <= 3) ) {
        Bm = get_ybgcorr1_hs(l, npt, sphr->dm, hrl, fr, sphr->rDm1,
            B2, vcr, 1, &Bc, &Bv, &fcorr);
        if (l >= 2) lambda[l - 2] = fcorr;
      } else {
        for ( i = 0; i < npt; i++ ) vcr[i] = vck[i] = 0;
        fcorr = 0;
        Bm = Bc;
      }
      sphr_r2k(sphr, hrl, hk[l]);
      /* additional corrections */
      for ( i = 0; i < npt; i++ ) {
        vck[i] *= fcorr;
        wkl[i] += vck[i];
        wk[l][i] = wkl[i];
        wr[l][i] += vcr[i];
        yrl[i] += vcr[i];
      }
      hk0[l] += integre(npt, vcr, fr, sphr->rDm1);
    }

    /* only for the hard-sphere model */
    yd[l] = (yrl[sphr->dm-1] + yrl[sphr->dm])/2;

    By = get_zerosep(wr[l], sphr->ri)*l/(l+1);

    savevir(fnvir, dim, l+2, Bc, Bv, Bm, 0, 0, By, B2, mkcorr, fcorr);
  }
  savevirtail(fnvir, clock() - t1);

  sphr_close(sphr);
  FREE1DARR(fr,   npt);
  FREE1DARR(fk,   npt);
  FREE1DARR(rdfr, npt);

  FREE1DARR(yrl, npt);
  FREE1DARR(hrl, npt);
  FREE1DARR(wkl, npt);
  FREE2DARR(wr, nmax - 1, npt);
  FREE2DARR(hk, nmax - 1, npt);
  FREE1DARR(yd, nmax - 1);
  FREE1DARR(hk0, nmax - 1);

  FREE1DARR(vcr, npt);
  FREE1DARR(vck, npt);
  FREE2DARR(wk, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE1DARR(lambda, nmax - 1);

  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, rmax, Rmax, ffttype);
  return 0;
}
