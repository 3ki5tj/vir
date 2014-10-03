/* Computing t_l(r) from c_l(r) computed from the OZ relation
 * adapted from ievir.c */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"



int dim = D;
int nmax = 12;
xdouble T = 1;
xdouble beta = 1;
double rmax = 0;
xdouble Rmax = 0;
int numpt = 32768;
int ffttype = 1;
int verbose = 0;
char *fncrtr = NULL;

int smoothpot = 0; /* smooth potential */
int gaussf = 0; /* Gaussian model */
int invexp = 0; /* inverse potential */
int dolj = 0; /* Lennard-Jones potential */
char systitle[32];

int doyr = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "computing the virial coefficients from the PY/HNC closure for hard-sphere fluids";
  argopt_add(ao, "-D", "%d", &dim, "dimension integer");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (exact)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type, 0: integral grid points, 1: half-integer grid points");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "--yr", "%b", &doyr, "deduce the bridge diagram from y(r)");
#ifdef DHT
  argopt_add(ao, "--disk", "%d", &dhtdisk, "use disk for discrete Hankel transform (DHT)");
  argopt_add(ao, "--tmpdir", "%s", &slowdht_tmpdir, "directory for temporary files");
  argopt_add(ao, "--dhtblock", "%u", &slowdht_block, "block size for DHT input/output");
  argopt_add(ao, "--dhtquiet", "%b", &slowdht_quiet, "suppress DHT output");
#endif
  argopt_add(ao, "-G", "%b", &gaussf, "Gaussian model instead of hard spheres");
  argopt_add(ao, "-I", "%d", &invexp, "exponent of the inverse potential r^(-n)");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);

  if ( dim < 1 ) argopt_help(ao);
  if ( dim % 2 == 0 ) {
#ifndef DHT
    fprintf(stderr, "cannot do even dimensions %d, define DHT\n", dim);
    exit(1);
#endif
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



/* load the reference c(r) that comes from Mayer sampling */
__inline double *loadcrraw(const char *fn, int *npt, double **ri)
{
  double *cr = NULL;
  int dim, order, i;
  char s[80];
  FILE *fp;

  *ri = NULL;

  xfopen(fp, fn, "r", return NULL);

  /* read the number of points */
  if ( fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "bad file %s\n", fn);
    goto EXIT;
  }
  if ( sscanf(s + 1, "%d%d%d", &dim, &order, npt) != 3 ) {
    fprintf(stderr, "cannot load information from %s\n%s", fn, s);
    goto EXIT;
  }

  /* load the array */
  xnew(cr, *npt);
  xnew(*ri, *npt);
  for ( i = 0; i < *npt; i++ ) {
    if ( fgets(s, sizeof s, fp) == NULL ) {
      fprintf(stderr, "%s corrupted for i %d\n", fn, i);
      free(*ri);
      free(cr);
      cr = NULL;
      goto EXIT;
    }
    sscanf(s, "%lf%lf", (*ri) + i, cr + i);
  }
EXIT:
  fclose(fp);
  return cr;
}



/* given the array (xi, yi), evaluate the value at x */
__inline double interp(double x, double *xi, double *yi,
    int imin, int imax, int cutimin, int cutimax)
{
  int i;
  double gam;

  for ( i = imin; i < imax; i++ ) if ( xi[i] > x ) break;
  if ( i == imax ) { /* extrapolate */
    if ( cutimax ) return 0;
    i--;
  } else if ( i == imin ) { /* extrapolate */
    if ( cutimin ) return 0;
    i++;
  }
  /* linear interpolation */
  gam = (xi[i] - x) / (xi[i] - xi[i-1]);
  return gam * yi[i-1] + (1 - gam) * yi[i];
}



#define loadcr(l, npt, crl, ri) loadcorl("cr", l, npt, crl, ri)

#define loadyrtr(l, npt, crl, ri) loadcorl("yrtr", l, npt, crl, ri)

/* load c(r) from Mayer sampling */
__inline int loadcorl(const char *prefix,
    int l, int npt, xdouble *crl, xdouble *ri)
{
  double *ri0 = NULL, *cr0 = NULL, r, y;
  int n0 = 0, dm0 = 0, imin0, imax0, cutimin0, cutimax0;
  char fn[80];
  int i;

  /* load the raw data */
  sprintf(fn, "%sD%dn%d.dat", prefix, dim, l+2);
  cr0 = loadcrraw(fn, &n0, &ri0);
  if ( cr0 == NULL )
    return -1;
  /* get the hard core boundary */
  for ( dm0 = 0; dm0 < n0; dm0++ )
    if ( ri0[dm0] > 1 ) break;

  for ( i = 0; i < npt; i++ ) {
    r = (double) ri[i];
    /* determine the interpolation range */
    if ( r < 1 ) {
      imin0 = 0;
      imax0 = dm0;
      cutimin0 = 1;
      cutimax0 = 0;
    } else {
      imin0 = dm0;
      imax0 = n0;
      cutimin0 = 0;
      cutimax0 = 1;
    }
    y = interp(r, ri0, cr0, imin0, imax0, cutimin0, cutimax0);
    crl[i] = y;
  }
  free(cr0);
  free(ri0);
  return 0;
}



/* compute w(r) from y(r) */
__inline static void getwr(int l, int npt, xdouble *wr, xdouble **yr)
{
  int i, u;
  xdouble *a, *b;

  xnew(a, l + 1);
  xnew(b, l + 1);
  for ( i = 0; i < npt; i++) {
    for ( a[0] = 1, u = 1; u <= l; u++ ) a[u] = yr[u][i];
    log_series(l + 1, a, b);
    wr[i] = b[l];
  }
  free(a);
  free(b);
}



static int intgeq(int nmax, int npt, xdouble rmax, xdouble Rmax, int ffttype)
{
  xdouble Bc, B2;
  xdouble *fr, *rdfr = NULL;
  xdouble *crl, *trl, *tkl;
  xdouble **yr = NULL, *yrtrl = NULL, *wrl = NULL;
  xdouble **ck;
  int i, l, l0 = 1;
  sphr_t *sphr;

  (void) Rmax;
  sphr = sphr_open(dim, npt, rmax, Rmax, ffttype);
  B2 = gaussf ? sphr->B2g : sphr->B2hs;
  /* TODO: compute B2 for inverse potential */
  if ( invexp > 0 ) B2 = 0;
  printf("D %d, dr %f, dm %d, rmax %f, ffttype %d, B2 %.10e\n",
      dim, (double) sphr->rmax/npt, sphr->dm, (double) sphr->rmax,
      ffttype, (double) B2);

  MAKE1DARR(fr, npt);
  if ( smoothpot ) MAKE1DARR(rdfr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE1DARR(yrtrl, npt);
  MAKE1DARR(wrl, npt);
  MAKE2DARR(yr, nmax - 1, npt);
  MAKE1DARR(tkl, npt);
  MAKE2DARR(ck, nmax - 1, npt);

  /* compute f(r) */
  mkfr(npt, beta, NULL, fr, rdfr, sphr->ri, sphr->dm, gaussf, invexp, dolj);
  COPY1DARR(crl, fr, npt)

  for ( i = 0; i < npt; i++ )
    yr[0][i] = 1;

  for ( l = l0; l < nmax - 1; l++ ) {
    if ( l == 2 ) {
      /* this is the PY solution, exact for l - 1 = 1 */
      for ( i = 0; i < npt; i++ )
        crl[i] = trl[i] * fr[i];
    }
    /* try to load c_{l - 1}(r), starts with c_2(r) */
    if ( l >= 3 && loadcr(l-1, npt, crl, sphr->ri) != 0 ) {
      fprintf(stderr, "stop at c_%d(r)\n", l - 1);
      break;
    }
    if (doyr) {
      /* try to load y_{l-1}(r), starts with y_2(r) */
      if ( l >= 3 && loadyrtr(l-1, npt, yrtrl, sphr->ri) != 0 ) {
        fprintf(stderr, "stop at y_%d(r)\n", l - 1);
        break;
      }
      for ( i = 0; i < npt; i++ )
        yr[l-1][i] = trl[i] + (l >= 3 ? yrtrl[i] : 0);
      getwr(l - 1, npt, wrl, yr);
      for ( i = 0; i < npt; i++ ) /* bridge function */
        yrtrl[i] = wrl[i] - trl[i];
    }
    Bc = get_Bc(l-1, npt, crl, sphr->rDm1);
    printf("n%4d, Bc %20.12" XDBLPRNF "e\n", l+1, Bc/pow_si(B2, l));
    savecrtr(fncrtr, l-1, npt, sphr->ri, crl, trl, yrtrl, wrl);

    /* c_l(r) --> c_l(k) for the previous l */
    sphr_r2k(sphr, crl, ck[l-1]);

    /* compute t_l(k) from c_0(k), ... c_{l-1}(k) */
    get_tk_oz(l, npt, ck, tkl);

    /* t_l(k) --> t_l(r) */
    sphr_k2r(sphr, tkl, trl);
  }

  sphr_close(sphr);
  FREE1DARR(fr,   npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(crl,  npt);
  FREE1DARR(trl,  npt);
  FREE1DARR(yrtrl, npt);
  FREE1DARR(wrl, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(tkl,  npt);
  FREE2DARR(ck, nmax - 1, npt);

  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, rmax, Rmax, ffttype);
  return 0;
}
