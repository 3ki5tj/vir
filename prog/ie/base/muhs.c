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
xdouble rho = (xdouble) 0.7L;
xdouble drho = (xdouble) 0.05L;
int itmax = 10000;
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
xdouble hcs = 1.83436;
char *fnBr; /* bridge function */



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing pressure and its derivatives from the PY closure";
  argopt_add(ao, "-r", "%" XDBLSCNF "f", &rho, "maximal density rho");
  argopt_add(ao, "--drho", "%" XDBLSCNF "f", &drho, "rho step size");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-d", "%" XDBLSCNF "f", &delta, "delta lambda");
  argopt_add(ao, "-O", "%b", &savecr, "save cr and tr");
  argopt_add(ao, "-s", "%" XDBLSCNF "f", &sqrs, "parameter of the quadratic closure");
  argopt_add(ao, "--py", "%b", &dopy, "do the PY closure");
  argopt_add(ao, "--hnc", "%b", &dohnc, "do the HNC closure");
  argopt_add(ao, "--ir", "%b", &doir, "inverse Rowlinson closure");
  argopt_add(ao, "--hc", "%b", &dohc, "Hutchinson-Conkie closure");
  argopt_add(ao, "--sc", "%b", &dosc, "do the self-consistent closure based on P");
  argopt_add(ao, "--SC", "%b", &doSC, "do the self-consistent closure based on mu");
  argopt_add(ao, "--Br", NULL, &fnBr, "file for the bridge function");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( dopy ) sqrs = 0;
  printf("rmax %f, rho %g, delta %g\n",
      (double) rmax, (double) rho, (double) delta);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
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
        fprintf(stderr, "%s: no line %d, l %d\n", fn, i, l);
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



/* combine the bridge function */
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
  if ( doir ) return irs += ds;
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
    sphr(npt, cr, ck, facr2k, plan, arr, ri, ki);
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
      cr[i] = x;
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
      dcr[i] = x;
    }
    //printf("round %d, errmax %g\n", it, errmax); getchar();
    if ( errmax < 1e-7 ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



/* compute the parameter for the self-consistent closure based on pressure */
static xdouble correct(int npt, xdouble rho, int dm,
    xdouble *cr, xdouble *tr, xdouble *fr, xdouble *ri2,
    xdouble *dcr, xdouble *dtr, xdouble surfr)
{
  int i;
  xdouble num = 0, den = 0, num2 = 0, den2 = 0, y, Dy, w, Dw, ds;

  (void) dcr;

  /* 1. compute s */
  for ( i = 0; i < npt; i++ ) {
    num += cr[i] * ri2[i];
    getyr(tr[i], NULL, &w, NULL);
    //w = 0.5 * tr[i] * tr[i];
    den += (1 + fr[i]) * w * ri2[i];
  }
  // checking code
  //y = 1+tr[dm]+.5*sqrs*tr[dm]*tr[dm];
  //Dy = 1 + sqrs*tr[dm];
  //w = 0.5*tr[dm]*tr[dm];
  //Dw = tr[dm];
  y = getyr(tr[dm], &Dy, &w, &Dw);
  num2 = (y + Dy * dtr[dm] * rho * .5) * surfr / dim;
  den2 = (w + Dw * dtr[dm] * rho * .5) * surfr / dim;
  num += num2;
  den += den2;
  ds = -num/den;
  //printf("num %g,%g, den %g,%g dm %d, y %g, %g, w %g, %g, s %g\n", (double) num, (double) num2, (double) den, (double) den2, dm,
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



/* solve the lambda case, lambda is the charging parameter */
static void iterc(int npt, xdouble rho, xdouble *ri, xdouble *ki,
    xdouble *Cr, xdouble *Tr, xdouble *Ck, xdouble *Tk,
    xdouble *Br, xdouble lambda,
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
        Yr = EXP(Tr[i] + Br[i]*lambda*lambda);
      } else {
        Yr = getyr(Tr[i], NULL, NULL, NULL);
      }
      x = (1 + Fr[i]) * Yr - (1 + Tr[i]);
      if ((err = FABS(Cr[i] - x)) > errmax) errmax = err;
      Cr[i] = x;
    }
    //printf("round %d, errmax %g, cr %g, %g\n", it, (double) errmax, (double) cr[0], (double) Cr[0]); getchar();
    if ( errmax < tol ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(lambda) functions, where lambda is the charging parameter */
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
      dcr[i] = x;
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
static xdouble getmu(int npt, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *br,
    xdouble *dcr, xdouble *dtr, xdouble *ri2,
    xdouble *mu1, xdouble *mu2,
    xdouble *ratio, xdouble *dmu)
{
  int i;
  xdouble yr, dyr, mu0, Bi, dBi, Bh = 0, hdB = 0, Bdh = 0;

  /* beta mu = Int { -c + B + (1/2) t h } + ... */
  mu0 = 0;
  *dmu = 0;
  for ( i = 0; i < npt; i++ ) {
    yr = getyr(tr[i], &dyr, NULL, NULL);
    Bi = LOG(yr) - tr[i];
    dBi = (dyr/yr - 1) * dtr[i];
    hdB += dBi*(cr[i] + tr[i]) * ri2[i];
    Bdh += Bi*(dcr[i] + dtr[i]) * ri2[i];

    if (br != 0) {
      Bi = dBi = br[i];
    }
    mu0 += (-cr[i] + Bi + .5*tr[i]*(cr[i] + tr[i])) * ri2[i];
    *dmu += (-dcr[i] + dBi + dtr[i]*tr[i] + .5*dtr[i]*cr[i]
            + .5*tr[i]*dcr[i] + dBi*(cr[i]+tr[i])) * ri2[i];
    //*dmu += (-dcr[i] + lambda * fr[i] * dtr[i]) * ri2[i];
    Bh += Bi*(cr[i] + tr[i]) * ri2[i];
  }
  if ( FABS(*ratio) < 1e-6 ) {
    *ratio = dohnc ? 0: hdB/(hdB + Bdh);
  }
  mu0 *= rho;
  Bh *= rho;
  *dmu *= rho;
  *mu1 = mu0 + Bh*2/3; /* low density limit */
  *mu2 = mu0 + Bh*(*ratio);
  return mu0;
}



static int integ(int npt, xdouble rmax, xdouble rhomax, xdouble rhodel)
{
  xdouble dr, dk, facr2k, fack2r, surfr, surfk, *bphi;
  xdouble mu0, mu0a, mu0b, rab0, dmu0, mu1, mu1a, mu1b, rab1, dmu1;
  xdouble *fr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk, *Dcr, *Dtr, *Dck, *Dtk;
  xdouble *Fr, *Cr, *Tr, *Ck, *Tk;
  xdouble *dCr, *dTr, *dCk, *dTk;
  xdouble *arr, *ri, *ki, *ri2, *ki2, s = 0, ds;
  int i, dm, sci;
  FFTWPFX(plan) plan = NULL;
  char fnout[80] = "iemu.dat";
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
  xnew(Cr, npt);
  xnew(Tr, npt);
  xnew(Ck, npt);
  xnew(Tk, npt);
  xnew(dCr, npt);
  xnew(dTr, npt);
  xnew(dCk, npt);
  xnew(dTk, npt);

  /* solvent-solvent interaction */
  for ( dm = -1, i = 0; i < npt; i++ ) {
    fr[i] = (ri[i] < 1) ? -1 : 0;
    Fr[i] = (ri[i] < 1) ? -(1-delta) : 0;
    if ( dm < 0 && ri[i] > 1) dm = i;
  }

  /* open the report file */
  sprintf(fnout, "iemuhs%s%s.dat",
      dohnc?"hnc":dopy?"py":dohc?"hc":doir?"ir":"sqr",
      dosc?"sc":doSC?"SC":"");
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

    /* 1. solve the case of lambda = 1 */
    iter(npt, rho, ri, ki, cr, tr, ck, tk, Br,
        fr, facr2k, fack2r, plan, arr, itmax);

    /* self-consistently determine s based on the pressure
     * it appears that this only makes the result slightly worse */
    if (dosc) {
      for (sci = 0; ; sci++) {
        /* differentiating with respect to rho */
        iterd(npt, rho, ri, ki, Dcr, Dtr, Dck, Dtk, tr, ck, tk,
              fr, facr2k, fack2r, plan, arr, itmax);
        ds = correct(npt, rho, dm, cr, tr, fr, ri2, Dcr, Dtr, surfr);
        iter(npt, rho, ri, ki, cr, tr, ck, tk, Br,
            fr, facr2k, fack2r, plan, arr, itmax);
        s = updates(ds);
        //printf("s %g, ds %g\n", (double) s, (double) ds); getchar();
        if (FABS(ds) < 1e-4) break;
      }
      //fprintf(stderr, "self-consistent iteration finished in %d rounds, s %g\n", sci, (double) s);
    }

    /* 2. compute dc(r), dt(r), differentiation w.r.t. lambda, at lambda = 1 */
    iterdc(npt, rho, ri, ki, dcr, dtr, dck, dtk, cr, tr, ck, tk,
        fr, fr, facr2k, fack2r, plan, arr, itmax);

    if (doSC) {
      for (sci = 0; ; sci++) {
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

    rab0 = 0;
    mu0 = getmu(npt, rho, cr, tr, Br, dcr, dtr, ri2, &mu0a, &mu0b, &rab0, &dmu0);

    if ( savecr )
      output(npt, ri, cr, tr, dcr, dtr, fr, bphi, "cr.dat");

    /* 3. solve the case of lambda = 1 - delta  */
    iterc(npt, rho, ri, ki, Cr, Tr, Ck, Tk, Br, 1-delta, cr, tr, ck, tk,
        Fr, facr2k, fack2r, plan, arr, itmax);

    /* 4. compute dc(r), dt(r), differentiation w.r.t. lambda, at lambda = 1 - delta */
    iterdc(npt, rho, ri, ki, dCr, dTr, dCk, dTk, Cr, Tr, Ck, Tk,
        Fr, fr, facr2k, fack2r, plan, arr, itmax);

    rab1 = rab0;
    mu1 = getmu(npt, rho, Cr, Tr, Br, dCr, dTr, ri2, &mu1a, &mu1b, &rab1, &dmu1);

    if ( savecr )
      output(npt, ri, Cr, Tr, dCr, dTr, Fr, bphi, "cr1.dat");

    printf("rho %5.3f, mu%9.4f, mua%9.4f, mub%9.4f, rab %6.4f "
           "dmu%9.4f,%9.4f,%9.4f;%9.4f; s %g\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b, (double) rab0,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s);
    fprintf(fp, "%8.6f %10.6f %10.6f %10.6f %10.8f %10.6f %10.6f %10.6f %10.6f %g\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b, (double) rab0,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s);
  }

  free(bphi);
  free(fr);
  free(cr);
  free(tr);
  free(ck);
  free(tk);
  free(dcr);
  free(dtr);
  free(dck);
  free(dtk);

  free(Fr);
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

