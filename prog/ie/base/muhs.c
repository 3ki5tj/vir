/* chemical potential related quantities
 * for the hard-sphere potential */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <fftw3.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"


#ifndef QUAD
#define LDBL /* use long double by default */
#endif

#include "xdouble.h"

#ifdef NOFFTW
#define XDOUBLE xdouble
#include "fft.h"
typedef void *FFTWPFX(plan);
#else
#include <fftw3.h>
#endif



int dim = 3; /* currently default to 3 */
int numpt = 8192;
xdouble rmax = (xdouble) 20.48L;
xdouble T = (xdouble) 10;
xdouble beta;
xdouble rho = (xdouble) 0.7L;
xdouble drho = (xdouble) 0.05L;
int itmax = 10000;
xdouble damp = 1;
xdouble tol = (xdouble) 1e-8L;
xdouble delta = (xdouble) 0.0001L;
int savecr = 0;
int verbose = 0;

int dopy = 0;
int dohnc = 0;
int doir = 0;
int dohc = 0;
int dosc = 0;
int doSC = 0;
xdouble sqrs = 0.16564;
xdouble irs = 0.16464;
xdouble hcs = 0.83436;

char *fnBr; /* bridge function */
xdouble slope = 0.1;

enum { SYS_HS, SYS_LJ };
int sys = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing chemical potential from integral equations";
  argopt_add(ao, "-r", "%" XDBLSCNF "f", &rho, "maximal density rho");
  argopt_add(ao, "--drho", "%" XDBLSCNF "f", &drho, "rho step size");
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &delta, "delta xi");
  argopt_add(ao, "-p", "%" XDBLSCNF "f", &damp, "damping factor for iteration");
  argopt_add(ao, "-O", "%b", &savecr, "save cr and tr");
  argopt_add(ao, "-s", "%" XDBLSCNF "f", &sqrs, "parameter of the quadratic closure");
  argopt_add(ao, "--py", "%b", &dopy, "do the PY closure");
  argopt_add(ao, "--hnc", "%b", &dohnc, "do the HNC closure");
  argopt_add(ao, "--ir", "%b", &doir, "inverse Rowlinson closure");
  argopt_add(ao, "--hc", "%b", &dohc, "Hutchinson-Conkie closure");
  argopt_add(ao, "--sc", "%b", &dosc, "do the self-consistent closure based on P");
  argopt_add(ao, "--SC", "%b", &doSC, "do the self-consistent closure based on mu");
  argopt_add(ao, "--Br", NULL, &fnBr, "file for the bridge function");
  argopt_add(ao, "-k", "%" XDBLSCNF "f", &slope, "k in gamma = 2/3 - k * rho");
  argopt_add(ao, "--lj", "%b", &sys, "do the Lennard-Jones fluid");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1./T;
  if ( dopy ) sqrs = 0;
  if ( doir ) irs = sqrs;
  if ( dohc ) hcs = sqrs;
  printf("T %g, rmax %f, rho %g, delta %g\n", (double) T,
      (double) rmax, (double) rho, (double) delta);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* return the potential phi(r), and -r*phi'(r)*/
static xdouble potlj(xdouble r, xdouble sig, xdouble eps, xdouble *ndphir)
{
  xdouble invr6 = (sig*sig)/(r*r), u;

  invr6 = invr6*invr6*invr6;
  u = 4*invr6*(invr6 - 1);
  if (u > 1000) u = 1000;
  *ndphir = invr6*(48*invr6 - 24);
  if (*ndphir > 12000) *ndphir = 12000;
  *ndphir *= eps;
  return eps*u;
}



/* given the array (xi, yi), evaluate the value at x */
__inline static double interp(double x, double *xi, double *yi,
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



/* load components of the bridge function */
__inline static int getBrs(const char *fn, int lmax, int npt,
    xdouble *ri, xdouble **Brs)
{
  FILE *fp;
  char buf[8192];
  int l, i, id = 0, imin, numpt;
  double *r, cr, tr, *br, wr;

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fn);
    return -1;
  }
  if (fgets(buf, sizeof buf, fp) == NULL ||
      buf[0] != '#') {
    fprintf(stderr, "%s: no information line\n", fn);
    return -1;
  }
  sscanf(buf + 1, "%d", &numpt);
  xnew(r, numpt);
  xnew(br, numpt);
  for ( l = 1; l <= lmax; l++ ) {
    /* read the raw data */
    for ( i = 0; i < numpt; i++ ) {
      if ( fgets(buf, sizeof buf, fp) == NULL ) {
        if (i != 0) {
          fprintf(stderr, "%s: no line %d, l %d\n", fn, i, l);
        } else {
          fprintf(stderr, "%s stop at order l %d\n", fn, l - 1);
        }
        goto EXIT;
      }
      if ( 6 != sscanf(buf, "%lf %lf %lf %d %lf %lf",
                       &r[i], &cr, &tr, &id, &br[i], &wr) ||
          id != l ) {
        fprintf(stderr, "%s: error on line %d, l %d\n%s", fn, i, l, buf);
        goto EXIT;
      }
    }
    fgets(buf, sizeof buf, fp);
    if ( !isspace(buf[0]) ) {
      fprintf(stderr, "error at the end of %s\n", fn);
      break;
    }

    if ( l == 1 ) continue;

    /* interpolation */
    for ( imin = 0; imin < numpt; imin++ )
      if (r[imin] > 0.1) break;
    for ( i = 0; i < npt; i++ )
      Brs[l][i] = interp((double) ri[i], r, br, imin, numpt, 0, 1);
  }
EXIT:
  free(r);
  free(br);
  fclose(fp);
  return l - 1;
}



/* combine the bridge function from the density components */
__inline static int combBr(int lmax, int npt, xdouble *Br, xdouble **Brs, xdouble rho)
{
  int i, l;
  xdouble rhopow;

  for ( i = 0; i < npt; i++ ) Br[i] = 0;
  rhopow = rho;
  for ( l = 2; l <= lmax; l++ ) {
    rhopow *= rho;
    for ( i = 0; i < npt; i++ ) {
      Br[i] += Brs[l][i] * rhopow;
    }
  }
  return 0;
}



/* get the value at zero separation
 * compute y(x) at x = 0 */
/* compute
 *    out(k) = 2*fac/k Int {from 0 to infinity} in(r) r sin(k r) dr */
static void sphr(int npt, xdouble *in, xdouble *out, xdouble fac,
    FFTWPFX(plan) p, xdouble *arr, xdouble *ri, xdouble *ki)
{
  int i;

  for ( i = 0; i < npt; i++ ) /* form in(x) * x */
    arr[i] = in[i] * ri[i];
#ifdef NOFFTW
  sint11(arr, npt);
#else
  FFTWPFX(execute)(p);
#endif
  for ( i = 0; i < npt; i++ ) /* form out(k) / k */
    out[i] = arr[i] * fac / ki[i];
}



/* compute the cavity function */
__inline static xdouble getyr(xdouble tr, xdouble *dy, xdouble *w, xdouble *dw)
{
  xdouble u, z;

  if ( dohnc ) { /* HNC */
    z = EXP(tr);
    if (dy) *dy = z;
    if (w) *w = z - 1 - tr;
    if (dw) *dw = z - 1;
    return z;
  } else if ( doir ) { /* inverse Rowlinson */
    z = EXP(tr);
    if (dy) *dy = 1 + irs * (z - 1);
    if (w) *w = z - 1 - tr;
    if (dw) *dw = z - 1;
    return 1 + tr + irs * (z - 1 - tr);
  } else if ( dohc ) { /* Hutchinson-Conkie */
    u = 1 + hcs * tr;
    z = POW(u, 1./hcs);
    if (dy) *dy = z/u;
    if (w) *w = (-LOG(u) + 1 - 1/u) * z / (hcs * hcs);
    if (dw) *dw = (-LOG(u) + (1 - hcs) * (1 - 1/u)) * z / u / (hcs * hcs);
    return z;
  } else { /* quadratic */
    if (dy) *dy = 1 + sqrs * tr;
    if (w) *w = .5 * tr * tr;
    if (dw) *dw = tr;
    return 1 + tr + sqrs * .5 * tr * tr;
  }
}



__inline static xdouble updates(xdouble ds)
{
  if ( dohnc || dopy ) return 0;
  else if ( doir ) return irs += ds;
  else if ( dohc ) return hcs += ds;
  else return sqrs += ds;
}



static void iter(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *Br, xdouble *fr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, err, errmax, yr;

  for ( it = 0; it < itmax; it++ ) {
    sphr(npt, cr, ck, facr2k, plan, arr, ri, ki); /* c(r) --> c(k) */
    for ( i = 0; i < npt; i++ ) {
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    }
    sphr(npt, tk, tr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      if ( Br != 0 ) {
        yr = EXP(tr[i] + Br[i]);
      } else {
        yr = getyr(tr[i], NULL, NULL, NULL);
      }
      x = (1 + fr[i]) * yr - (1 + tr[i]);
      if ((err = FABS(cr[i] - x)) > errmax) errmax = err;
      cr[i] += damp * (x - cr[i]);
    }
    if ( errmax < tol ) break;
  }
  //printf("iter %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(rho) functions */
static void iterd(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    xdouble *fr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, hk, dy, err, errmax;

  for ( i = 0; i < npt; i++ ) dcr[i] = fr[i];
  for ( it = 0; it < itmax; it++ ) {
    sphr(npt, dcr, dck, facr2k, plan, arr, ri, ki);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dtk[i] = hk*hk + rho*hk*(2 + rho*hk)*dck[i];
    }
    sphr(npt, dtk, dtr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* dc = (1 + f) dy - dt = [(1 + f) Y' - 1] dt */
      getyr(tr[i], &dy, NULL, NULL);
      x = ((1 + fr[i]) * dy - 1) * dtr[i];
      if ((err = FABS(dcr[i] - x)) > errmax) errmax = err;
      dcr[i] += damp * (x - dcr[i]);
    }
    //printf("round %d, errmax %g\n", it, errmax); getchar();
    if ( errmax < 1e-7 ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



/* compute the parameter for the self-consistent closure based on pressure */
static xdouble correct(int npt, xdouble rho, int dm,
    xdouble *cr, xdouble *tr, xdouble *fr, xdouble *ri2,
    xdouble *dcr, xdouble *dtr,
    xdouble surfr, xdouble *rdfr)
{
  int i;
  xdouble num = 0, den = 0, y, Dy, w, Dw, ds;

  (void) dcr;

  /* 1. compute s */
  for ( i = 0; i < npt; i++ ) {
    /* Dy = dy / dt
     * w = dy / ds
     * Dw = d^2 y / ds dt */
    y = getyr(tr[i], &Dy, &w, &Dw);
    //w = 0.5 * tr[i] * tr[i];
    num += cr[i] * ri2[i];
    den += (1 + fr[i]) * w * ri2[i]; /* d_xi c = (1 + f) d_xi y */
    if ( sys != SYS_HS ) {
      num += rdfr[i] * (y + .5 * rho * Dy * dtr[i]) / dim;
      den += rdfr[i] * (w + .5 * rho * Dw * dtr[i]) / dim;
    }
  }
  if ( sys == SYS_HS ) {
    // checking code
    //y = 1+tr[dm]+.5*sqrs*tr[dm]*tr[dm];
    //Dy = 1 + sqrs*tr[dm];
    //w = 0.5*tr[dm]*tr[dm];
    //Dw = tr[dm];
    y = getyr(tr[dm], &Dy, &w, &Dw);
    num += (y + Dy * dtr[dm] * rho * .5) * surfr / dim;
    den += (w + Dw * dtr[dm] * rho * .5) * surfr / dim;
  }
  ds = -num/den;
  //printf("num %g, den %g, dm %d, y %g, %g, w %g, %g, s %g\n", (double) num, (double) den, dm,
  //    (double) y, (double) Dy, (double) w, (double) Dw, (double) sqrs); getchar();

  /* 2. use s to correct the correlation functions */
  for ( i = 0; i < npt; i++ ) {
    getyr(tr[i], NULL, &w, NULL);
    w *= dtr[i];
    cr[i] += (1 + fr[i]) * w * ds;
  }
  return ds;
}



/* compute the parameter for the self-consistent closure based on mu */
static xdouble correct2(int npt, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *fr, xdouble *ri2,
    xdouble *dcr, xdouble *dtr, xdouble *Dcr, xdouble *Dtr)
{
  int i;
  xdouble num = 0, den = 0, y, Dy, Dw, ds;

  (void) Dcr;

  /* 1. compute s */
  for ( i = 0; i < npt; i++ ) {
    y = getyr(tr[i], &Dy, NULL, &Dw);
    num += (dcr[i] - fr[i] * y - rho * fr[i] * Dy * Dtr[i]) * ri2[i];
    den += ((1 + fr[i]) * dtr[i] - rho * fr[i] * Dtr[i]) * Dw * ri2[i];
  }
  //printf("num %g, den %g, rho %g\n", (double) num, (double) den, (double) rho); getchar();
  ds = -num/den;

  /* 2. use s to correct the correlation functions */
  for ( i = 0; i < npt; i++ ) {
    getyr(tr[i], NULL, &Dy, NULL);
    cr[i] += (1 + fr[i]) * Dy * ds;
  }
  return ds;
}



/* solve the xi != 1 case, xi is the charging parameter */
static void iterc(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *Cr, xdouble *Tr, xdouble *Ck, xdouble *Tk,
    xdouble *Br, xdouble xi,
    const xdouble *cr, const xdouble *tr,
    const xdouble *ck, const xdouble *tk,
    xdouble *Fr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, err, errmax, Yr;

  (void) tr;
  for ( i = 0; i < npt; i++ ) Cr[i] = cr[i];
  for ( it = 0; it < itmax; it++ ) {
    /* get C(k) from C(r) */
    sphr(npt, Cr, Ck, facr2k, plan, arr, ri, ki);
    /* T(k) = rho C(k) h(k) */
    for ( i = 0; i < npt; i++ )
      Tk[i] = rho*Ck[i]*(ck[i] + tk[i]);
    /* get T(r) from T(k) */
    sphr(npt, Tk, Tr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* C(r) = (1 + F(r)) Y(r) - T(r) - 1 */
      if ( Br != NULL ) {
        Yr = EXP(Tr[i] + Br[i]*xi*xi);
      } else {
        Yr = getyr(Tr[i], NULL, NULL, NULL);
      }
      x = (1 + Fr[i]) * Yr - (1 + Tr[i]);
      if ((err = FABS(Cr[i] - x)) > errmax) errmax = err;
      Cr[i] += damp * (x - Cr[i]);
    }
    //printf("round %d, errmax %g, cr %g, %g\n", it, (double) errmax, (double) cr[0], (double) Cr[0]); getchar();
    if ( errmax < tol ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(xi) functions, where xi is the charging parameter */
static void iterdc(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *cr, const xdouble *tr, const xdouble *ck, const xdouble *tk,
    xdouble *fr, xdouble *dfr, xdouble facr2k, xdouble fack2r,
    FFTWPFX(plan) plan, xdouble *arr, int itmax)
{
  int i, it;
  xdouble x, err, errmax, yr, dyr;

  for ( i = 0; i < npt; i++ ) dcr[i] = cr[i];
  for ( it = 0; it < itmax; it++ ) {
    /* get dc(k) from dc(r) */
    sphr(npt, dcr, dck, facr2k, plan, arr, ri, ki);
    /* dt(k) = rho dc(k) h(k) */
    for ( i = 0; i < npt; i++ )
      dtk[i] = rho*dck[i]*(ck[i] + tk[i]);
    /* get dt(r) from dt(k) */
    sphr(npt, dtk, dtr, fack2r, plan, arr, ki, ri);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* dc(r) = -dt(r) + (1 + f(r)) dy(r) + df(r) y(r)
       * for the hard sphere we assume that df(r) = f(r) */
      yr = getyr(tr[i], &dyr, NULL, NULL);
      dyr *= dtr[i];
      x = -dtr[i] + (1 + fr[i]) * dyr + dfr[i] * yr;
      if ((err = FABS(dcr[i] - x)) > errmax) errmax = err;
      dcr[i] += damp * (x - dcr[i]);
    }
    //printf("round %d, errmax %g\n", it, (double) errmax); getchar();
    if ( errmax < tol ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



static void output(int npt, xdouble *ri,
    xdouble *cr, xdouble *tr, xdouble *dcr, xdouble *dtr,
    xdouble *fr, xdouble *bphi, char *fn)
{
  int i;
  FILE *fp;

  xfopen(fp, fn, "w", return);
  for ( i = 0; i < npt; i++ ) {
    fprintf(fp, "%8.6f %14.8f %14.8f %14.8f %14.8f %14.8f %g\n",
        (double) ri[i],
        (double) cr[i], (double) tr[i],
        (double) dcr[i], (double) dtr[i],
        (double) fr[i], (double) bphi[i]);
  }
  fclose(fp);
}



/* compute the excess chemical potential and its derivatives */
static xdouble getmu(int npt, xdouble rho, xdouble *fr,
    xdouble *cr, xdouble *tr, xdouble *br,
    xdouble *dcr, xdouble *dtr,
    xdouble *Dcr, xdouble *Dtr,
    xdouble *ffr, xdouble *ri2,
    xdouble *mu1, xdouble *mu2, xdouble *mu2r, xdouble *dmu,
    xdouble *mu3, xdouble *mu4, xdouble *mu5, xdouble *mu5th,
    xdouble *mu6, xdouble *mu6r1, xdouble *mu6r2)
{
  int i;
  xdouble yr, dyr, mu0, B, dBdt, dB;
  xdouble sf = 0, sc = 0, stt = 0, sth = 0, stth = 0, sfff = 0, sB = 0, sBh = 0;
  xdouble stdt = 0, scdt = 0, stdc = 0, sdB = 0, shdB = 0, sBdh = 0;
  xdouble sDc = 0, stDc = 0, scDt = 0, stDt = 0, sDB = 0, shDB = 0, sBDh = 0;
  xdouble corr1 = 0, dcorr1 = 0, Dcorr1 = 0, corr2 = 0, dcorr2 = 0, Dcorr2 = 0;
  xdouble numd, dend1, dend2, numD, denD1, denD2, Q0, Q1, r = 0, r1 = 0, r2 = 0, thc, thv;
  int ctype = dohnc ? 1 : 0;

  mu0 = 0;
  *dmu = 0;
  for ( i = 0; i < npt; i++ ) {
    yr = getyr(tr[i], &dyr, NULL, NULL);
    B = LOG(yr) - tr[i];
    dBdt = dyr/yr - 1;
    dB = dBdt * dtr[i];

    sB += B * ri2[i];
    sBh += B * (cr[i] + tr[i]) * ri2[i];
    sdB += dB * ri2[i];
    shdB += dB * (cr[i] + tr[i]) * ri2[i];
    sBdh += B * (dcr[i] + dtr[i]) * ri2[i];

    sc += cr[i] * ri2[i];
    sf += fr[i] * ri2[i];
    sth += tr[i] * (cr[i] + tr[i]) * ri2[i];
    stt += tr[i] * tr[i] * ri2[i];
    stth += tr[i] * tr[i] * (cr[i] + tr[i]) * ri2[i];
    scdt += cr[i] * dtr[i] * ri2[i]; /* scdt should be equal to stdc */
    stdc += tr[i] * dcr[i] * ri2[i];
    stdt += tr[i] * dtr[i] * ri2[i];
    if ( ffr != NULL ) sfff += ffr[i] * fr[i] * ri2[i];

    if (br != 0) { /* explicit bridge function */
      B = br[i];
      dB = 2*B;
    }
    /* beta dmu = Int ( -dc + dB + (1/2) dt h + (1/2) t dh + dB h) dr */
    *dmu += rho * (-dcr[i] + dB + dtr[i]*tr[i] + .5*dtr[i]*cr[i]
            + .5*tr[i]*dcr[i] + dB*(cr[i]+tr[i])) * ri2[i];

    if ( ctype == 1 )
      dcorr2 += tr[i] * (tr[i] * dcr[i] + (2 * cr[i] + 3 * tr[i]) * dtr[i]) * ri2[i];

    if ( Dcr != NULL && Dtr != NULL ) {
      sDc += Dcr[i] * ri2[i];
      stDc += tr[i] * Dcr[i] * ri2[i];
      scDt += cr[i] * Dtr[i] * ri2[i];
      stDt += tr[i] * Dtr[i] * ri2[i];
      sBDh += B * (Dcr[i] + Dtr[i]) * ri2[i];
      shDB += dBdt * Dtr[i] * (cr[i] + tr[i]) * ri2[i];
      sDB += dBdt * Dtr[i] * ri2[i];
      if ( ctype == 1 )
        Dcorr2 += tr[i] * (tr[i] * Dcr[i] + (2 * cr[i] + 3 * tr[i]) * Dtr[i]) * ri2[i];
    }
  }
  /* beta mu = Int ( -c + B + (1/2) t h ) dr
   *         + Int {0 to 1} dxi Int dB h dr */
  mu0 = rho * (-sc + sB + .5 * sth);

  if ( ctype == 0 ) {
    corr1 = sB;
    dcorr1 = sdB;
    Dcorr1 = sDB;
    corr2 = sBh;
    dcorr2 = shdB + sBdh;
    Dcorr2 = shDB + sBDh;
  } else if ( ctype == 1 ) {
    corr1 = stt;
    dcorr1 = 2*stdt;
    Dcorr1 = 2*stDt;
    corr2 = stth;
  }

  numd = shdB; /* residue d/dxi */
  dend1 = dcorr1;
  dend2 = dcorr2;
  if ( FABS(*mu2r) < 1e-6 ) { /* if mu2r is not give */
    *mu2r = dohnc ? 0: numd/dcorr2;
  }
  *mu1 = rho * (-sc + sB + .5 * sth + corr2 * 2./3); /* low density limit */
  *mu2 = rho * (-sc + sB + .5 * sth + corr2 * (*mu2r));
  if ( Dcr != NULL && Dtr != NULL ) {
    /* -Dcr = -gr D(tr+Br) + Dtr = -(hr Dtr + hr DBr + DBr)
     *      = -hk Dtk - hr DBr - DBr
     * (1/2) D(hr tr) = (1/2) D(hk tk) = (1/2) (Dhk tk + Dtk hk)
     * -Dcr + (1/2) D (hr tr)
     *  = (1/2) (tk Dhk - hk Dtk) - hr DBr - DBr
     *  = (1/2) (tk Dck - ck Dtk) - hr DBr - DBr
     *  = (1/2) (tr Dcr - cr Dtr) - hr DBr - DBr
     * So,
     * D { rho Int [-cr + (1/2) tr hr] }
     * = Int [ -cr + (1/2) tr hr ]
     *   + rho Int [ (1/2) (tr Dcr - cr Dtr) - hr DBr - DBr]
     * Or
     * D { rho Int [-cr + Br + (1/2) tr hr] }
     * = Int [ -cr + Br + (1/2) tr hr]
     *   + Int [ (1/2) (tr Dcr - cr Dtr) - hr DBr]
     * */
    numD = -(sB + 0.5*sth) + rho*(0.5*(scDt - stDc) + shDB); /* residue d/drho */
    denD1 = corr1 + rho * Dcorr1;
    denD2 = corr2 + rho * Dcorr2;
    thc = numD / denD1;
    /* local power expansion */
    /*
    *mu3 = 0;
    for ( i = 0; i < npt; i++ ) {
      if ( ffr != NULL ) {
        Q0 = cr[i] - fr[i] - rho * ffr[i] * fr[i];
        Q1 = Dcr[i] - ffr[i] * fr[i];
        *mu3 += -rho * ri2[i] * (fr[i] + fr[i]*ffr[i]*rho/2 + Q0*Q0/(rho*Q1 + Q0 + 1e-8));
      } else {
        Q0 = cr[i] - fr[i];
        Q1 = Dcr[i];
        *mu3 += -rho * ri2[i] * (fr[i] + Q0 * Q0/ (Q0 + Q1*rho + 1e-8));
      }
    }
    */
    /* global/local polynomial approximation */
    if ( ffr != NULL )
      *mu3 = -rho * ((sf + sc)/2 + rho*(sfff - sDc)/12);
    else
      *mu3 = -rho * (2*sc/3 + sf/3 - rho*sDc/6);

    /* global power approximation */
    Q0 = sc - sf - rho * sfff;
    Q1 = sDc - sfff;
    *mu4 = -rho * (sf + sfff*rho/2 +  Q0 * Q0 / (rho*Q1 + Q0 + 1e-8));
  } else {
    *mu4 = *mu3 = -rho * (sc + sf)/2;
  }
  if (mu5th != NULL) *mu5th = thc;
  *mu5 = rho * (-sc + sB + 0.5 * sth + thc * corr1);
  r1 = (numd * denD2 - numD * dend2) / (dend1 * denD2 - denD1 * dend2 + 1e-8);
  r2 = (numd * denD1 - numD * dend1) / (dend2 * denD1 - denD2 * dend1 - 1e-8);
  if (mu6r1 != NULL) *mu6r1 = r1;
  if (mu6r2 != NULL) *mu6r2 = r2;
  *mu6 = rho * (-sc + sB + 0.5 * sth + r1 * corr1 + r2 * corr2);
  return mu0;
}



/* compute the pressure */
static xdouble getpres(int npt, xdouble rho, xdouble *fr,
    xdouble *cr, xdouble *tr,
    xdouble *Dcr, xdouble *Dtr, xdouble *ri2,
    xdouble *ck, xdouble *ki2,
    xdouble *pres1, xdouble *pres2, xdouble *p2r,
    xdouble *pres3, xdouble *pres4, xdouble *ffr)
{
  int i;
  xdouble spk = 0, spr = 0, x, yr, dyr, w, dwdt, Dw;
  xdouble sc = 0, stc = 0, sDc = 0, scDt = 0, stDc = 0, sf = 0, sfff = 0;
  xdouble sew = 0, setw = 0, seDw = 0, sewDt = 0, setDw = 0, r = 0, corr = 0, Dcorr = 0;
  int ctype = 1;

  /* k-space sum */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if (FABS(x) < 1e-8) {
      spk += -x*x*x*(1./3 + x*.25) * ki2[i];
    } else {
      spk += (log(1 - x) + x + x*x*.5) * ki2[i];
    }
  }

  /* real-space sum */
  for ( i = 0; i < npt; i++ ) {
    sf += fr[i] * ri2[i];
    sfff += ffr[i] * fr[i] * ri2[i];
    sc += cr[i] * ri2[i];
    stc += cr[i] * tr[i] * ri2[i];
    yr = getyr(tr[i], &dyr, NULL, NULL);
    w = yr - tr[i] - 1;
    dwdt = dyr - 1;
    sew += (fr[i] + 1) * w * ri2[i];
    setw += (fr[i] + 1) * tr[i] * w * ri2[i];
    if ( Dcr != NULL && Dtr != NULL ) {
      sDc += Dcr[i] * ri2[i];
      scDt += cr[i] * Dtr[i] * ri2[i];
      stDc += tr[i] * Dcr[i] * ri2[i];
      Dw = dwdt * Dtr[i];
      sewDt += (fr[i] + 1) * w * Dtr[i] * ri2[i];
      seDw += (fr[i] + 1) * Dw * ri2[i];
      setDw += (fr[i] + 1) * tr[i] * Dw * ri2[i];
    }
  }
  spr = -0.5 * rho * rho * (sc - stc);
  if ( ctype == 0 ) {
    corr = .5 * rho * rho * sew;
    Dcorr = rho * sew + .5 * rho * rho * seDw;
  } else {
    corr = .5 * rho * rho * (sew + setw);
    Dcorr = rho * (sew + setw) + .5 * rho * rho * (seDw + sewDt + setDw);
  }

  *pres1 = rho + spk + spr + 0.5 * corr;
  if ( Dcr != NULL && Dtr != NULL ) {
    xdouble num = .5 * rho * rho * (sDc - scDt + stDc);
    xdouble den = Dcorr;
    xdouble Q0, Q1;
    if (!dopy) r = num/den;
    //printf("%g/%g r %g\n", (double)num, (double)den, (double)r);
    Q0 = sc - sf - rho * sfff;
    Q1 = sDc - sfff;
    *pres3 = rho - rho*rho*(sf*.5 + sfff*rho/3 +  Q0*Q0/(rho*Q1 + 2*Q0 + 1e-8));
    *pres4 = rho - rho*rho*(sf*0.15 + sc*0.35 + rho*(sfff/30 - sDc/20));
  }
  if (p2r != NULL) *p2r = r;
  *pres2 = rho + spk + spr + r * corr;
  return rho + spr + spk;
}



static int integ(int npt, xdouble rmax, xdouble rhomax, xdouble rhodel)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk, *bphi;
  xdouble mu0, mu0a, mu0b, rab0, dmu0, mu1, mu1a, mu1b, rab1, dmu1;
  xdouble mu0c, mu1c, mucth = 0, mu0d, mu1d, mueth = 0, mu0e, mu1e, mu0f, mu1f;
  xdouble mu0fr1 = 0, mu0fr2 = 0, mu1fr1 = 0, mu1fr2 = 0;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk, *Dcr, *Dtr, *Dck, *Dtk;
  xdouble *Fr, *rdFr, *Cr, *Tr, *Ck, *Tk;
  xdouble *ffr, *Ffr;
  xdouble *dCr, *dTr, *dCk, *dTk;
  xdouble pres0, pres0a, pres0b, rpres0 = 0, pres0c, pres0d;
  xdouble *arr, *ri, *ki, *ri2, *ki2, s = sqrs, ds;
  int i, dm, sci;
  FFTWPFX(plan) plan = NULL;
  char fnout[80] = "iemu.dat", buf[80] = "";
  FILE *fp;
  int Brlmax = -1;
  xdouble *Br = NULL, **Brs = NULL;

  xnew(arr, npt);
  xnew(ri, npt);
  xnew(ki, npt);
  xnew(ri2, npt);
  xnew(ki2, npt);
  plan = FFTWPFX(plan_r2r_1d)(npt, arr, arr, FFTW_RODFT11, FFTW_ESTIMATE);

  dr = rmax / npt;
  dk = PI / (dr * npt);

  facr2k = PI*2 * dr;
  fack2r = pow_si(PI*2, -2) * dk;
  for ( i = 0; i < npt; i++ ) {
    ri[i] = dr * (i * 2 + 1) / 2;
    ki[i] = dk * (i * 2 + 1) / 2;
  }

  surfr = PI*4;
  surfk = surfr * pow_si(PI*2, -3);
  for ( i = 0; i < npt; i++ ) {
    ri2[i] = surfr * ri[i] * ri[i] * dr;
    ki2[i] = surfk * ki[i] * ki[i] * dk;
  }

  xnew(bphi, npt);
  xnew(fr, npt);
  xnew(rdfr, npt);
  xnew(cr, npt);
  xnew(tr, npt);
  xnew(ck, npt);
  xnew(tk, npt);
  xnew(dcr, npt);
  xnew(dtr, npt);
  xnew(dck, npt);
  xnew(dtk, npt);
  xnew(Dcr, npt);
  xnew(Dtr, npt);
  xnew(Dck, npt);
  xnew(Dtk, npt);

  xnew(Fr, npt);
  xnew(rdFr, npt);
  xnew(Cr, npt);
  xnew(Tr, npt);
  xnew(Ck, npt);
  xnew(Tk, npt);
  xnew(dCr, npt);
  xnew(dTr, npt);
  xnew(dCk, npt);
  xnew(dTk, npt);

  xnew(ffr, npt);
  xnew(Ffr, npt);

  /* solvent-solvent interaction */
  for ( dm = -1, i = 0; i < npt; i++ ) {
    if ( sys == SYS_HS ) {
      fr[i] = (ri[i] < 1) ? -1 : 0;
    } else {
      xdouble x;
      bphi[i] = beta * potlj(ri[i], 1, 1, &rdfr[i]);
      x = EXP(-bphi[i]);
      fr[i] = x - 1;
      rdfr[i] *= beta * x;
    }
    Fr[i] = fr[i] * (1 - delta);
    rdFr[i] = rdfr[i] * (1 - delta);
    if ( dm < 0 && ri[i] > 1) dm = i;
  }

  /* compute the fff */
  sphr(npt, fr, tk, facr2k, plan, arr, ri, ki); /* f(r) --> t(k) */
  for ( i = 0; i < npt; i++ ) tk[i] = tk[i] * tk[i];
  sphr(npt, tk, ffr, fack2r, plan, arr, ki, ri); /* t(k) --> t(r) */
  sphr(npt, Fr, Tk, facr2k, plan, arr, ri, ki); /* F(r) --> T(k) */
  for ( i = 0; i < npt; i++ ) Tk[i] = Tk[i] * Tk[i];
  sphr(npt, Tk, Ffr, fack2r, plan, arr, ki, ri); /* T(k) --> T(r) */

  /* open the report file */
  if ( sys != SYS_HS ) sprintf(buf, "T%g", (double) T);
  else buf[0] = '\0';
  sprintf(fnout, "iemu%s%s%s%s.dat",
      sys?"lj":"hs", buf,
      fnBr?"Br":dohnc?"hnc":dopy?"py":dohc?"hc":doir?"ir":"sqr",
      dosc?"scp":doSC?"scmu":"");
  if ((fp = fopen(fnout, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fnout);
    return -1;
  }
  fprintf(fp, "# %g %g\n", (double) rho, (double) drho);

  if ( fnBr != NULL ) {
    xnew(Br, npt);
    Brlmax = 7;
    xnew(Brs, Brlmax);
    xnew(Brs[0], Brlmax * npt);
    for ( i = 0; i < 2*npt; i++ ) Brs[0][i] = 0;
    for ( i = 1; i < Brlmax; i++ ) Brs[i] = Brs[0] + npt*i;
    Brlmax = getBrs(fnBr, Brlmax, npt, ri, Brs);
    if (Brlmax > 0) {
      fprintf(stderr, "loaded the bridge function from %s\n", fnBr);
    } else {
      fprintf(stderr, "failed to load the bridge function from %s\n", fnBr);
      exit(1);
    }
  }

  /* initialize c(r) for iteration */
  for ( i = 0; i < npt; i++ ) cr[i] = fr[i];

  for ( rho = rhodel; rho < rhomax + tol; rho += rhodel ) {
    if (Br != NULL) combBr(Brlmax, npt, Br, Brs, rho);

    /* 1. solve the case of xi = 1 */
    iter(npt, rho, ri, ki, cr, tr, ck, tk, Br,
        fr, facr2k, fack2r, plan, arr, itmax);

    /* 2. differentiating with respect to rho */
    iterd(npt, rho, ri, ki, Dcr, Dtr, Dck, Dtk, tr, ck, tk,
          fr, facr2k, fack2r, plan, arr, itmax);

    /* self-consistently determine s based on the pressure
     * it appears that this only makes the result slightly worse */
    if (dosc) {
      for (sci = 1; ; sci++) {
        /* differentiating with respect to rho */
        iterd(npt, rho, ri, ki, Dcr, Dtr, Dck, Dtk, tr, ck, tk,
              fr, facr2k, fack2r, plan, arr, itmax);
        ds = correct(npt, rho, dm, cr, tr, fr, ri2, Dcr, Dtr, surfr, rdfr);
        iter(npt, rho, ri, ki, cr, tr, ck, tk, Br,
            fr, facr2k, fack2r, plan, arr, itmax);
        s = updates(ds);
        //printf("s %g, ds %g\n", (double) s, (double) ds); getchar();
        if (FABS(ds) < 1e-4) break;
      }
      if (verbose)
        fprintf(stderr, "self-consistent iteration finished in %d rounds, s %g\n", sci, (double) s);
    }

    /* 3. compute dc(r), dt(r), differentiation w.r.t. xi, at xi = 1 */
    iterdc(npt, rho, ri, ki, dcr, dtr, dck, dtk, cr, tr, ck, tk,
        fr, fr, facr2k, fack2r, plan, arr, itmax);

    if (doSC) {
      for (sci = 1; ; sci++) {
        /* differentiating with respect to rho */
        iterd(npt, rho, ri, ki, Dcr, Dtr, Dck, Dtk, tr, ck, tk,
              fr, facr2k, fack2r, plan, arr, itmax);
        ds = correct2(npt, rho, cr, tr, fr, ri2, dcr, dtr, Dcr, Dtr);
        iter(npt, rho, ri, ki, cr, tr, ck, tk, Br,
            fr, facr2k, fack2r, plan, arr, itmax);
        iterdc(npt, rho, ri, ki, dcr, dtr, dck, dtk, cr, tr, ck, tk,
            fr, fr, facr2k, fack2r, plan, arr, itmax);
        s = updates(ds);
        if (FABS(ds) < 1e-4) break;
      }
      //fprintf(stderr, "self-consistent iteration finished in %d rounds, s %g\n", sci, (double) s);
    }

    rab0 = (Br != NULL) ? (2./3 - slope * rho) : 0;
    mu0 = getmu(npt, rho, fr, cr, tr, Br, dcr, dtr, Dcr, Dtr, ffr, ri2,
        &mu0a, &mu0b, &rab0, &dmu0,
        &mu0c, &mu0d, &mu0e, &mueth, &mu0f, &mu0fr1, &mu0fr2);
    pres0 = getpres(npt, rho, fr, cr, tr, Dcr, Dtr, ri2, ck, ki2,
        &pres0a, &pres0b, &rpres0, &pres0c, &pres0d, ffr);

    if ( savecr )
      output(npt, ri, cr, tr, dcr, dtr, fr, bphi, "cr.dat");

    /* 4. solve the case of xi = 1 - delta  */
    iterc(npt, rho, ri, ki, Cr, Tr, Ck, Tk, Br, 1-delta, cr, tr, ck, tk,
        Fr, facr2k, fack2r, plan, arr, itmax);

    /* 5. compute dc(r), dt(r), differentiation w.r.t. xi, at xi = 1 - delta */
    iterdc(npt, rho, ri, ki, dCr, dTr, dCk, dTk, Cr, Tr, Ck, Tk,
        Fr, fr, facr2k, fack2r, plan, arr, itmax);

    rab1 = rab0;
    mu1 = getmu(npt, rho, Fr, Cr, Tr, Br, dCr, dTr, NULL, NULL, Ffr, ri2,
        &mu1a, &mu1b, &rab1, &dmu1,
        &mu1c, &mu1d, &mu1e, NULL, &mu1f, &mu1fr1, &mu1fr2);

    if ( savecr )
      output(npt, ri, Cr, Tr, dCr, dTr, Fr, bphi, "cr1.dat");

    printf("rho %5.3f, muv%9.4f,%9.4f,%9.4f rab %6.4f "
           "dmu%9.4f,%9.4f,%9.4f;%9.4f s%8.5f\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b,
           (double) rab0,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s);
    printf("rho %5.3f, muc%9.4f,%9.4f,%9.4f th%+8.4f "
           "mu %9.4f,%9.4f,%9.4f\n",
           (double) rho,
           (double) mu0c, (double) mu0d, (double) mu0e, (double) mueth,
           (double) mu0f, (double) mu0fr1, (double) mu0fr2);
    printf("rho %5.3f, P  %9.4f,%9.4f,%9.4f rP %7.4f Pcd%9.4f,%9.4f\n",
           (double) rho,
           (double) pres0, (double) pres0a, (double) pres0b, (double) rpres0,
           (double) pres0c, (double) pres0d);
    fprintf(fp, "%8.6f "
                "%10.6f %10.6f %10.6f %10.8f "
                "%12.6f %12.6f %12.6f %12.6f "
                "%10.8f "
                "%10.6f %10.8f %10.6f %10.6f "
                "%10.6f %10.8f %10.8f "
                "%10.6f %10.6f %10.6f %10.8f %10.6f %10.6f\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b, (double) rab0,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s,
           (double) mu0c, (double) mu0d, (double) mu0e, (double) mueth,
           (double) mu0f, (double) mu0fr1, (double) mu0fr2,
           (double) pres0, (double) pres0a, (double) pres0b, (double) rpres0,
           (double) pres0c, (double) pres0d);
  }

  free(bphi);
  free(fr);
  free(rdfr);
  free(cr);
  free(tr);
  free(ck);
  free(tk);
  free(dcr);
  free(dtr);
  free(dck);
  free(dtk);

  free(Fr);
  free(rdFr);
  free(Cr);
  free(Tr);
  free(Ck);
  free(Tk);
  free(dCr);
  free(dTr);
  free(dCk);
  free(dTk);

  free(arr);
  free(ri);
  free(ki);
  free(ri2);

  free(Br);
  if (Brs != NULL) {
    free(Brs[0]);
    free(Brs);
  }

  fclose(fp);
  fprintf(stderr, "saved report to %s\n", fnout);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);

  integ(numpt, rmax, rho, drho);
#ifndef NOFFTW
  FFTWPFX(cleanup)();
#endif
  return 0;
}

