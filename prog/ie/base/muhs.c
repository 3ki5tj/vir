/* chemical potential related quantities
 * for the hard-sphere potential */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"


/* warning: routines in this file may not be stable */
//#include "nlfit.h"



int dim = D; /* currently default to 3 */
int numpt = 8192;
int ffttype = 1;
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

int dohs = 1;
int gaussf = 0;
int invexp = 0;
int dolj = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing chemical potential from integral equations";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "--rho", "%" XDBLSCNF "f", &rho, "maximal density rho");
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
  argopt_add(ao, "-G", "%b", &gaussf, "do the Gaussian fluid");
  argopt_add(ao, "--invexp", "%d", &invexp, "exponent e of the r^{-e} fluid");
  argopt_add(ao, "--lj", "%b", &dolj, "do the Lennard-Jones fluid");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  beta = 1./T;
  if ( dopy ) sqrs = 0;
  if ( doir ) irs = sqrs;
  if ( dohc ) hcs = sqrs;
  if ( gaussf || invexp || dolj ) dohs = 0;
  printf("T %g, rmax %f, rho %g, delta %g\n", (double) T,
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
    if ( fgets(buf, sizeof buf, fp) == NULL
      || !isspace(buf[0]) ) {
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
__inline static int combBr(int lmax, int npt, xdouble *Br, xdouble *dBr,
    xdouble **Brs, xdouble rho)
{
  int i, l;
  xdouble rhopow;

  for ( i = 0; i < npt; i++ ){
    dBr[i] = Br[i] = 0;
    rhopow = rho;
    for ( l = 2; l <= lmax; l++ ) {
      dBr[i] += Brs[l][i] * l * rhopow;
      rhopow *= rho;
      Br[i] += Brs[l][i] * rhopow;
    }
  }
  return 0;
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



/* solve the integral equation */
static void iter(sphr_t *sphr, xdouble rho,
    xdouble *cr, xdouble *tr, xdouble *ck, xdouble *tk,
    xdouble *Br, xdouble *fr, int itmax)
{
  int i, npt = sphr->npt, it;
  xdouble x, err, errmax, yr;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, cr, ck); /* c(r) --> c(k) */
    for ( i = 0; i < npt; i++ ) {
      tk[i] = rho * ck[i] * ck[i] / (1 - rho * ck[i]);
    }
    sphr_k2r(sphr, tk, tr);
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
static void iterd(sphr_t *sphr, xdouble rho,
    xdouble *dcr, xdouble *dtr, xdouble *dck, xdouble *dtk,
    const xdouble *tr, const xdouble *ck, const xdouble *tk,
    const xdouble *Br, const xdouble *dBr, const xdouble *fr,
    int itmax)
{
  int i, npt = sphr->npt, it;
  xdouble x, hk, dydt, err, errmax;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, dcr, dck);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dtk[i] = hk*hk + rho*hk*(2 + rho*hk)*dck[i];
    }
    sphr_r2k(sphr, dtk, dtr);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* dc = (1 + f) dy - dt = [(1 + f) Y' - 1] dt */
      if (Br) {
        x = (1 + fr[i]) * EXP(tr[i] + Br[i]) * (dtr[i] + dBr[i]) - dtr[i];
      } else {
        getyr(tr[i], &dydt, NULL, NULL);
        x = ((1 + fr[i]) * dydt - 1) * dtr[i];
      }
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
    xdouble B2, xdouble *rdfr)
{
  int i;
  xdouble num = 0, den = 0, y, dydt, w, dwdt, ds;

  (void) dcr;

  /* 1. compute s */
  for ( i = 0; i < npt; i++ ) {
    /* dydt = dy / dt
     * w = dy / ds
     * dwdt = d^2 y / ds dt */
    y = getyr(tr[i], &dydt, &w, &dwdt);
    //w = 0.5 * tr[i] * tr[i];
    num += cr[i] * ri2[i];
    den += (1 + fr[i]) * w * ri2[i]; /* d_xi c = (1 + f) d_xi y */
    if ( !dohs ) {
      num += rdfr[i] * (y + .5 * rho * dydt * dtr[i]) / dim;
      den += rdfr[i] * (w + .5 * rho * dwdt * dtr[i]) / dim;
    }
  }
  if ( dohs ) {
    // checking code
    //y = 1+tr[dm]+.5*sqrs*tr[dm]*tr[dm];
    //dydt = 1 + sqrs*tr[dm];
    //w = 0.5*tr[dm]*tr[dm];
    //dwdt = tr[dm];
    y = getyr(tr[dm], &dydt, &w, &dwdt);
    num += (y + dydt * dtr[dm] * rho * .5) * 2 * B2;
    den += (w + dwdt * dtr[dm] * rho * .5) * 2 * B2;
  }
  ds = -num/den;
  //printf("num %g, den %g, dm %d, y %g, %g, w %g, %g, s %g\n", (double) num, (double) den, dm,
  //    (double) y, (double) dydt, (double) w, (double) dwdt, (double) sqrs); getchar();

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
  xdouble num = 0, den = 0, y, dydt, w, dwdt, ds;

  (void) Dcr;

  /* 1. compute s */
  for ( i = 0; i < npt; i++ ) {
    y = getyr(tr[i], &dydt, NULL, &dwdt);
    num += (dcr[i] - fr[i] * y - rho * fr[i] * dydt * Dtr[i]) * ri2[i];
    den += ((1 + fr[i]) * dtr[i] - rho * fr[i] * Dtr[i]) * dwdt * ri2[i];
  }
  //printf("num %g, den %g, rho %g\n", (double) num, (double) den, (double) rho); getchar();
  ds = -num/den;

  /* 2. use s to correct the correlation functions */
  for ( i = 0; i < npt; i++ ) {
    getyr(tr[i], NULL, &w, NULL);
    cr[i] += (1 + fr[i]) * w * ds;
  }
  return ds;
}



/* solve the xi != 1 case, xi is the charging parameter */
static void iterc(sphr_t *sphr, xdouble rho,
    xdouble *Cr, xdouble *Tr, xdouble *Ck, xdouble *Tk,
    xdouble *Br, xdouble xi,
    const xdouble *ck, const xdouble *tk,
    xdouble *Fr, int itmax)
{
  int i, npt = sphr->npt, it;
  xdouble x, err, errmax, Yr;

  for ( it = 0; it < itmax; it++ ) {
    /* get C(k) from C(r) */
    sphr_r2k(sphr, Cr, Ck);
    /* T(k) = rho C(k) h(k) */
    for ( i = 0; i < npt; i++ )
      Tk[i] = rho*Ck[i]*(ck[i] + tk[i]);
    /* get T(r) from T(k) */
    sphr_k2r(sphr, Tk, Tr);
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



/* solve d/d(rho) functions with xi != 1 */
static void itercd(sphr_t *sphr, xdouble rho,
    xdouble *dCr, xdouble *dTr, xdouble *dCk, xdouble *dTk,
    const xdouble *dck, const xdouble *dtk,
    const xdouble *Ck, const xdouble *Tr,
    const xdouble *ck, const xdouble *tk,
    xdouble *Fr, int itmax)
{
  int i, npt = sphr->npt, it;
  xdouble x, hk, dhk, dY, err, errmax;

  for ( it = 0; it < itmax; it++ ) {
    sphr_r2k(sphr, dCr, dCk);
    for ( i = 0; i < npt; i++ ) {
      hk = ck[i] + tk[i];
      dhk = dck[i] + dtk[i];
      dTk[i] = Ck[i]*hk + rho*hk*dCk[i] + rho*Ck[i]*dhk;
    }
    sphr_k2r(sphr, dTk, dTr);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* dC = (1 + F) dY - dT = [(1 + F) Y' - 1] dT */
      getyr(Tr[i], &dY, NULL, NULL);
      x = ((1 + Fr[i]) * dY - 1) * dTr[i];
      if ((err = FABS(dCr[i] - x)) > errmax) errmax = err;
      dCr[i] += damp * (x - dCr[i]);
    }
    //printf("round %d, errmax %g\n", it, errmax); getchar();
    if ( errmax < 1e-7 ) break;
  }
  //printf("it %d errmax %g\n", it, (double) errmax);
}



/* solve d/d(xi) functions, where xi is the charging parameter
 * which can be 1 or not 1 */
static void itercD(sphr_t *sphr, xdouble rho,
    xdouble *Dcr, xdouble *Dtr, xdouble *Dck, xdouble *Dtk,
    const xdouble *cr, const xdouble *tr, const xdouble *ck, const xdouble *tk,
    xdouble *fr, xdouble *Dfr, int itmax)
{
  int i, npt = sphr->npt, it;
  xdouble x, err, errmax, yr, Dyr;

  (void) cr;
  for ( it = 0; it < itmax; it++ ) {
    /* Dc(r) --> Dc(k) */
    sphr_r2k(sphr, Dcr, Dck);
    /* Dt(k) = rho Dc(k) h(k) */
    for ( i = 0; i < npt; i++ )
      Dtk[i] = rho*Dck[i]*(ck[i] + tk[i]);
    /* Dt(k) --> Dt(r) */
    sphr_k2r(sphr, Dtk, Dtr);
    for ( errmax = 0, i = 0; i < npt; i++ ) {
      /* Dc(r) = -Dt(r) + (1 + f(r)) Dy(r) + Df(r) y(r)
       * for the hard sphere we assume that Df(r) = f(r) */
      yr = getyr(tr[i], &Dyr, NULL, NULL);
      Dyr *= Dtr[i];
      x = -Dtr[i] + (1 + fr[i]) * Dyr + Dfr[i] * yr;
      if ((err = FABS(Dcr[i] - x)) > errmax) errmax = err;
      Dcr[i] += damp * (x - Dcr[i]);
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
    xdouble *cr, xdouble *tr,
    xdouble *Dcr, xdouble *Dtr, /* derivative w.r.t. xi */
    xdouble *dcr, xdouble *dtr, /* derivative w.r.t. rho */
    xdouble *br, xdouble *dbr,
    xdouble *ffr, xdouble *ri2,
    xdouble *mu1, xdouble *mu2, xdouble *mu2r, xdouble *dmu,
    xdouble *mu3, xdouble *mu4, xdouble *mu5, xdouble *mu5th,
    xdouble *mu6, xdouble *mu6r1, xdouble *mu6r2)
{
  int i;
  xdouble yr, dyr, mu0, B, h, dBdt, DB, Dh, dB, dh;
  xdouble sf = 0, sc = 0, stc = 0, stt = 0, sth = 0, stth = 0, sfff = 0, sB = 0, sBh = 0;
  xdouble sDc = 0, stDt = 0, scDt = 0, stDc = 0, sDB = 0, shDB = 0, sBDh = 0;
  xdouble sdc = 0, stdc = 0, scdt = 0, stdt = 0, sdB = 0, shdB = 0, sBdh = 0;
  xdouble corr1 = 0, Dcorr1 = 0, dcorr1 = 0, corr2 = 0, Dcorr2 = 0, dcorr2 = 0;
  xdouble numD, denD1, denD2, numd = 0, dend1 = 0, dend2 = 0, det;
  int ctype = dohnc ? 1 : 0;

  mu0 = 0;
  *dmu = 0;
  for ( i = 0; i < npt; i++ ) {
    if (br != 0) { /* explicit bridge function */
      B = br[i];
      DB = 2*B;
    } else {
      yr = getyr(tr[i], &dyr, NULL, NULL);
      dBdt = dyr/yr - 1;
      B = LOG(yr) - tr[i];
      DB = dBdt * Dtr[i];
    }
    h = cr[i] + tr[i];
    Dh = Dcr[i] + Dtr[i];

    sB += B * ri2[i];
    sBh += B * h * ri2[i];
    sDB += DB * ri2[i];
    shDB += h * DB * ri2[i];
    sBDh += B * Dh * ri2[i];

    sc += cr[i] * ri2[i];
    sf += fr[i] * ri2[i];
    stc += tr[i] * cr[i] * ri2[i];
    sth += tr[i] * h * ri2[i];
    stt += tr[i] * tr[i] * ri2[i];
    stth += tr[i] * tr[i] * h * ri2[i];
    sDc += Dcr[i] * ri2[i];
    scDt += cr[i] * Dtr[i] * ri2[i]; /* scDt should be equal to stDc */
    stDc += tr[i] * Dcr[i] * ri2[i];
    stDt += tr[i] * Dtr[i] * ri2[i];
    if ( ffr != NULL ) {
      sfff += ffr[i] * fr[i] * ri2[i];
    }

    if ( ctype == 1 ) {
      //corr2 += -.5 * tr[i] * tr[i] * h * ri2[i];
      //Dcorr2 += -.5 * tr[i] * (tr[i]*Dh + 2*h*Dtr[i]) * ri2[i];
      //corr2 += -.5 * tr[i] * tr[i]*fr[i] * ri2[i];
      //Dcorr2 += -.5 * tr[i] * (tr[i]*fr[i] + 2*fr[i]*Dtr[i]) * ri2[i];
      corr2 += -.5 * rho * tr[i]*ffr[i]*fr[i] * ri2[i];
      Dcorr2 += -.5 * rho * (2*tr[i] + Dtr[i])*ffr[i]*fr[i] * ri2[i];
      //corr2 += -.5*rho*rho * ffr[i]*ffr[i]*fr[i] * ri2[i];
      //Dcorr2 += -.5*rho*rho * (3*ffr[i]*ffr[i]*fr[i]) * ri2[i];
    }

    /* derivatives w.r.t. density */
    if ( dcr != NULL && dtr != NULL ) {
      if (br != 0) {
        dB = dbr[i];
      } else {
        dB = dBdt * dtr[i];
      }
      dh = dcr[i] + dtr[i];
      sdc += dcr[i] * ri2[i];
      stdc += tr[i] * dcr[i] * ri2[i];
      scdt += cr[i] * dtr[i] * ri2[i];
      stdt += tr[i] * dtr[i] * ri2[i];
      sdB += dB * ri2[i];
      sBdh += B * dh * ri2[i];
      shdB += dB * h * ri2[i];
      if ( ctype == 1 ) {
        //dcorr2 += -.5*tr[i]*(tr[i]*dh + 2*h*dtr[i]) * ri2[i];
        //dcorr2 += -.5*tr[i]*(2*fr[i]*dtr[i]) * ri2[i];
        dcorr2 += -.5*(tr[i] + rho*dtr[i])*ffr[i]*fr[i] * ri2[i];
        //dcorr2 += -.5*(2*rho*ffr[i]*ffr[i]*fr[i]) * ri2[i];
      }
    }
  }
  /* beta mu = Int ( -c + B + (1/2) t h ) dr
   *         + Int {0 to 1} Di Int DB h dr */
  mu0 = rho * (-sc + sB + .5 * sth);
  /* Differentation with respect to the charging parameter
   * beta dmu = Int ( -Dc + DB + h Dt h + h DB ) dr */
  *dmu = rho * (-sDc + sDB + scDt + stDt + shDB);

  if ( ctype == 0 ) {
    corr1 = sB;
    Dcorr1 = sDB;
    dcorr1 = sdB;
    corr2 = sBh;
    Dcorr2 = shDB + sBDh;
    //printf("hDB %g, BDh %g\n", (double)(shDB/sBh),(double)(sBDh/sBh));
    dcorr2 = shdB + sBdh;
    //printf("hdB %g, Bdh %g\n", (double)(rho*shdB/sBh),(double)(rho*sBdh/sBh));
  } else if ( ctype == 1 ) {
    corr1 = sc - (sf + rho * sfff);
    Dcorr1 = sDc - (sf + 2 * rho * sfff); /* d/dxi */
    dcorr1 = sdc - sfff; /* d/drho */
    //corr2 = sc;
    //Dcorr2 = sDc;
    //dcorr2 = sdc;
    //corr2 = stt;
    //Dcorr2 = 2*stDt;
    //dcorr2 = 2*stdt;
    //corr2 = stth;
    //corr1 = stt + corr2;
    //Dcorr1 = 2*stDt + Dcorr2;
    //dcorr1 = 2*stdt + dcorr2;
  }

  numD = shDB; /* residue d/Di */
  denD1 = Dcorr1;
  denD2 = Dcorr2;
  if ( FABS(*mu2r) < 1e-6 ) { /* if mu2r is not give */
    *mu2r = dohnc ? 0: numD/Dcorr2;
  }
  *mu1 = mu0 + rho * corr2 * (dohnc ? 0 : 2./3); /* low density limit */
  *mu2 = mu0 + rho * corr2 * (*mu2r); /* virial-route consistent */

  /* derivatives w.r.t. density */
  if ( dcr != NULL && dtr != NULL ) {
    xdouble Q0, Q1;

    /* -dcr = -gr d(tr+Br) + dtr = -(hr dtr + hr dBr + dBr)
     *      = -hk dtk - hr dBr - dBr
     * (1/2) d(hr tr) = (1/2) d(hk tk) = (1/2) (dhk tk + dtk hk)
     * -dcr + (1/2) d (hr tr)
     *  = (1/2) (tk dhk - hk dtk) - hr dBr - dBr
     *  = (1/2) (tk dck - ck dtk) - hr dBr - dBr
     *  = (1/2) (tr dcr - cr dtr) - hr dBr - dBr
     * So,
     * d { rho Int [-cr + (1/2) tr hr] }
     * = Int [ -cr + (1/2) tr hr ]
     *   + rho Int [ (1/2) (tr dcr - cr dtr) - hr dBr - dBr]
     * and
     * d { rho Int [-cr + Br + (1/2) tr hr] }
     * = Int [ -cr + Br + (1/2) tr hr]
     *   + Int [ (1/2) (tr dcr - cr dtr) - hr dBr]
     * */
    if ( FABS(*mu5th) < 1e-6 ) {
      numd = -(sB + 0.5*sth) + rho*(0.5*(scdt - stdc) + shdB); /* residue d/drho */
      dend1 = corr1 + rho * dcorr1;
      dend2 = corr2 + rho * dcorr2;
      *mu5th = numd / dend1;
    }
    /* local power expansion */
    /*
    *mu3 = 0;
    for ( i = 0; i < npt; i++ ) {
      xdouble Q0, Q1;
      if ( ffr != NULL ) {
        Q0 = cr[i] - fr[i] - rho * ffr[i] * fr[i];
        Q1 = dcr[i] - ffr[i] * fr[i];
        *mu3 += -rho * ri2[i] * (fr[i] + fr[i]*ffr[i]*rho/2 + Q0*Q0/(rho*Q1 + Q0 + 1e-8));
      } else {
        Q0 = cr[i] - fr[i];
        Q1 = dcr[i];
        *mu3 += -rho * ri2[i] * (fr[i] + Q0 * Q0/ (Q0 + Q1*rho + 1e-8));
      }
    }
    */
    /* global/local polynomial approximation */
    if ( ffr != NULL )
      *mu3 = -rho * ((sf + sc)/2 + rho*(sfff - sdc)/12);
    else
      *mu3 = -rho * (2*sc/3 + sf/3 - rho*sdc/6);

    /* global power approximation */
    Q0 = sc - sf - rho * sfff;
    Q1 = sdc - sfff;
    *mu4 = -rho * (sf + sfff*rho/2 +  Q0 * Q0 / (rho*Q1 + Q0 + 1e-8));
    //*mu4 = getmuiter(rho, sf, sfff, sc, sdc);
    //printf("rho %g, %g, %g\n", (double) rho, (double) *mu4, (double) mu4b); getchar();
  } else {
    *mu4 = *mu3 = -rho * (sc + sf)/2;
  }

  /* compressibility-route consistent */
  *mu5 = mu0 + rho * (*mu5th) * corr1;

  /* consistent for both the compressibility and the virial routes
   *    r1 * denD1 + r2 * denD2 = numD
   *    r1 * dend1 + r2 * dend2 = numd
   * */
  if ( FABS(*mu6r1) < 1e-6 && FABS(*mu6r2) < 1e-6 ) { /* if not give */
    det = denD1 * dend2 - dend1 * denD2 + 1e-8;
    *mu6r1 = (numD * dend2 - numd * denD2) / det;
    *mu6r2 = (numd * denD1 - numD * dend1) / det;
  }
  *mu6 = mu0 + rho * (*mu6r1 * corr1 + *mu6r2 * corr2);
  return mu0;
}



/* compute the pressure */
static xdouble getpres(int npt, xdouble rho, int dm, xdouble B2,
    xdouble *fr, xdouble *cr, xdouble *tr,
    xdouble *dcr, xdouble *dtr, xdouble *ri2,
    xdouble *ck, xdouble *ki2,
    xdouble *pres1, xdouble *pres2, xdouble *p2r,
    xdouble *pres3, xdouble *pres4, xdouble *ffr)
{
  int i;
  xdouble spk = 0, spr = 0, x, yr, dyr, w, dwdt, dw;
  xdouble sc = 0, stc = 0, sdc = 0, scdt = 0, stdc = 0, sf = 0, sfff = 0;
  xdouble sew = 0, setw = 0, sedw = 0, sewdt = 0, setdw = 0, r = 0, corr = 0, dcorr = 0;
  int ctype = 1;

  (void) dm;
  (void) B2;
  /* k-space sum */
  for ( i = 0; i < npt; i++ ) {
    x = rho * ck[i];
    if (FABS(x) < 1e-8) {
      spk += -x*x*x*(1./3 + x*.25) * ki2[i];
    } else {
      spk += (LOG(1 - x) + x + x*x*.5) * ki2[i];
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
    if ( dcr != NULL && dtr != NULL ) {
      sdc += dcr[i] * ri2[i];
      scdt += cr[i] * dtr[i] * ri2[i];
      stdc += tr[i] * dcr[i] * ri2[i];
      dw = dwdt * dtr[i];
      sewdt += (fr[i] + 1) * w * dtr[i] * ri2[i];
      sedw += (fr[i] + 1) * dw * ri2[i];
      setdw += (fr[i] + 1) * tr[i] * dw * ri2[i];
    }
  }
  spr = -0.5 * rho * rho * (sc - stc);
  if ( ctype == 0 ) {
    corr = .5 * rho * rho * sew;
    dcorr = rho * sew + .5 * rho * rho * sedw;
  } else {
    corr = .5 * rho * rho * (sew + setw);
    dcorr = rho * (sew + setw) + .5 * rho * rho * (sedw + sewdt + setdw);
  }

  *pres1 = rho + spk + spr + 0.5 * corr;
  if ( dcr != NULL && dtr != NULL ) {
    xdouble num = .5 * rho * rho * (sdc - scdt + stdc);
    xdouble den = dcorr;
    if (!dopy) r = num/den;
    *pres3 = rho - rho*rho*(sf*0.15 + sc*0.35 + rho*(sfff/30 - sdc/20));
    //printf("%g/%g r %g\n", (double)num, (double)den, (double)r);
    xdouble Q0, Q1;
    Q0 = sc - sf - rho * sfff;
    Q1 = sdc - sfff;
    *pres4 = rho - rho*rho*(sf*.5 + sfff*rho/3 +  Q0*Q0/(rho*Q1 + 2*Q0 + 1e-8));
    //*pres4 = getpresiter(rho, sf, sfff, sc, sdc);
  } else {
    *pres3 = *pres4 = 0;
  }
  if (p2r != NULL) *p2r = r;
  *pres2 = rho + spk + spr + r * corr;
  return rho + spr + spk;
  //return rho - rho*rho*(sc + B2*(1+tr[dm])*(1+tr[dm]));
}



static int integ(int npt, xdouble rmax, xdouble rhomax, xdouble rhodel)
{
  xdouble *bphi;
  xdouble mu0, mu0a, mu0b, mu0br, dmu0, mu1, mu1a, mu1b, dmu1;
  xdouble mu0c, mu1c, mu0d, mu1d, mu0e, mu1e, mu0f, mu1f;
  xdouble mueth, mufr1, mufr2;
  xdouble *fr, *rdfr, *cr, *tr, *ck, *tk;
  xdouble *dcr, *dtr, *dck, *dtk;
  xdouble *Dcr, *Dtr, *Dck, *Dtk;
  xdouble *Fr, *rdFr, *Cr, *Tr, *Ck, *Tk;
  xdouble *ffr, *Ffr;
  xdouble *dCr, *dTr, *dCk, *dTk;
  xdouble *DCr, *DTr, *DCk, *DTk;
  xdouble pres0, pres0a, pres0b, rpres0 = 0, pres0c, pres0d;
  xdouble s = sqrs, ds;
  int i, sci;
  char fnout[80] = "iemu.dat", buf[80] = "";
  FILE *fp;
  int Brlmax = -1;
  xdouble *Br = NULL, *dBr = NULL, **Brs = NULL;
  sphr_t *sphr;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);

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
  xnew(DCr, npt);
  xnew(DTr, npt);
  xnew(DCk, npt);
  xnew(DTk, npt);

  xnew(ffr, npt);
  xnew(Ffr, npt);

  mkfr(npt, beta, bphi, fr, rdfr, sphr->ri, sphr->dm, gaussf, invexp, dolj);
  for ( i = 0; i < npt; i++ ) {
    Fr[i] = fr[i] * (1 - delta);
    rdFr[i] = rdfr[i] * (1 - delta);
  }

  /* compute the fff */
  sphr_r2k(sphr, fr, tk); /* f(r) --> t(k) */
  for ( i = 0; i < npt; i++ ) tk[i] = tk[i] * tk[i];
  sphr_k2r(sphr, tk, ffr); /* t(k) --> t(r) */
  sphr_r2k(sphr, Fr, Tk); /* F(r) --> T(k) */
  for ( i = 0; i < npt; i++ ) Tk[i] = Tk[i] * Tk[i];
  sphr_k2r(sphr, Tk, Ffr); /* T(k) --> T(r) */

  /* open the report file */
  if ( !dohs ) sprintf(buf, "T%g", (double) T);
  else buf[0] = '\0';
  sprintf(fnout, "iemu%s%s%s%s.dat",
      dolj?"lj":gaussf?"gauss":invexp?"invexp":"hs", buf,
      fnBr?"Br":dohnc?"hnc":dopy?"py":dohc?"hc":doir?"ir":"sqr",
      dosc?"scp":doSC?"scmu":"");
  if ((fp = fopen(fnout, "w")) == NULL) {
    fprintf(stderr, "cannot open %s\n", fnout);
    return -1;
  }
  fprintf(fp, "# %g %g\n", (double) rho, (double) drho);

  if ( fnBr != NULL ) {
    xnew(Br, npt);
    xnew(dBr, npt);
    Brlmax = 7;
    xnew(Brs, Brlmax);
    xnew(Brs[0], Brlmax * npt);
    for ( i = 0; i < 2*npt; i++ ) Brs[0][i] = 0;
    for ( i = 1; i < Brlmax; i++ ) Brs[i] = Brs[0] + npt*i;
    Brlmax = getBrs(fnBr, Brlmax, npt, sphr->ri, Brs);
    if (Brlmax > 0) {
      fprintf(stderr, "loaded the bridge function from %s\n", fnBr);
    } else {
      fprintf(stderr, "failed to load the bridge function from %s\n", fnBr);
      exit(1);
    }
  }

  /* initialize c(r) for iteration */
  for ( i = 0; i < npt; i++ ) {
    Cr[i] = DCr[i] = dCr[i] = cr[i] = Dcr[i] = dcr[i] = fr[i];
  }

  for ( rho = rhodel; rho < rhomax + tol; rho += rhodel ) {
    if (Br != NULL) combBr(Brlmax, npt, Br, dBr, Brs, rho);

    /* 1. solve the case of xi = 1 */
    iter(sphr, rho, cr, tr, ck, tk, Br, fr, itmax);

    /* 2. differentiating with respect to rho */
    iterd(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, Br, dBr, fr, itmax);

    /* self-consistently determine s based on the pressure
     * it appears that this only makes the result slightly worse */
    if (dosc) {
      for (sci = 1; ; sci++) {
        /* differentiating with respect to rho */
        iterd(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, Br, dBr, fr, itmax);
        ds = correct(npt, rho, sphr->dm, cr, tr, fr, sphr->rDm1, dcr, dtr, sphr->B2hs, rdfr);
        iter(sphr, rho, cr, tr, ck, tk, Br, fr, itmax);
        s = updates(ds);
        //printf("s %g, ds %g\n", (double) s, (double) ds); getchar();
        if (FABS(ds) < 1e-4) break;
      }
      if (verbose)
        fprintf(stderr, "self-consistent iteration finished in %d rounds, s %g\n", sci, (double) s);
    }

    /* 3. compute Dc(r), Dt(r), differentiation w.r.t. xi, at xi = 1 */
    itercD(sphr, rho, Dcr, Dtr, Dck, Dtk, cr, tr, ck, tk, fr, fr, itmax);

    if (doSC) {
      for (sci = 1; ; sci++) {
        /* differentiating with respect to rho */
        iterd(sphr, rho, dcr, dtr, dck, dtk, tr, ck, tk, Br, dBr, fr, itmax);
        ds = correct2(npt, rho, cr, tr, fr, sphr->rDm1, Dcr, Dtr, dcr, dtr);
        iter(sphr, rho, cr, tr, ck, tk, Br, fr, itmax);
        itercD(sphr, rho, Dcr, Dtr, Dck, Dtk, cr, tr, ck, tk, fr, fr, itmax);
        s = updates(ds);
        if (FABS(ds) < 1e-4) break;
      }
      //fprintf(stderr, "self-consistent iteration finished in %d rounds, s %g\n", sci, (double) s);
    }

    /* compute thermodynamic quantities */
    mu0br = (Br != NULL) ? (2./3 - slope * rho) : 0;
    mufr1 = mufr2 = mueth = 0;
    mu0 = getmu(npt, rho, fr, cr, tr, Dcr, Dtr, dcr, dtr, Br, dBr, ffr, sphr->rDm1,
        &mu0a, &mu0b, &mu0br, &dmu0,
        &mu0c, &mu0d, &mu0e, &mueth, &mu0f, &mufr1, &mufr2);
    pres0 = getpres(npt, rho, sphr->dm, sphr->B2hs, fr, cr, tr, dcr, dtr,
        sphr->rDm1, ck, sphr->kDm1,
        &pres0a, &pres0b, &rpres0, &pres0c, &pres0d, ffr);

    if ( savecr )
      output(npt, sphr->ri, cr, tr, Dcr, Dtr, fr, bphi, "cr.dat");

    /* 4. solve the case of xi = 1 - delta  */
    iterc(sphr, rho, Cr, Tr, Ck, Tk, Br, 1-delta, ck, tk, Fr, itmax);

    /* 5. solve derivatives w.r.t. density, rho */
    itercd(sphr, rho, dCr, dTr, dCk, dTk, dck, dtk, Ck, Tr, ck, tk, Fr, itmax);

    /* 6. compute DC(r), DT(r), derivatives w.r.t. xi, at xi = 1 - delta */
    itercD(sphr, rho, DCr, DTr, DCk, DTk, Cr, Tr, Ck, Tk, Fr, fr, itmax);

    mu1 = getmu(npt, rho, Fr, Cr, Tr, DCr, DTr, dCr, dTr, Br, dBr, Ffr, sphr->rDm1,
        &mu1a, &mu1b, &mu0br, &dmu1,
        &mu1c, &mu1d, &mu1e, &mueth, &mu1f, &mufr1, &mufr2);

    if ( savecr )
      output(npt, sphr->ri, Cr, Tr, DCr, DTr, Fr, bphi, "Cr.dat");

    printf("rho %5.3f, muv%9.4f,%9.4f,%9.4f rab %6.4f "
           "dmu%9.4f,%9.4f,%9.4f;%9.4f s%8.5f\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b,
           (double) mu0br,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s);
    printf("rho %5.3f, muc%9.4f,%9.4f,%9.4f th%+8.4f "
           "mu %9.4f;%9.4f,%9.4f %9.4f\n",
           (double) rho,
           (double) mu0c, (double) mu0d, (double) mu0e, (double) mueth,
           (double) mu0f, (double) mufr1, (double) mufr2,
           (double) ((mu0f - mu1f)/delta));
    printf("rho %5.3f, P  %9.4f,%9.4f,%9.4f rP %7.4f Pcd%9.4f,%9.4f\n",
           (double) rho,
           (double) pres0, (double) pres0a, (double) pres0b, (double) rpres0,
           (double) pres0c, (double) pres0d);
    fprintf(fp, "%8.6f "
                "%10.6f %10.6f %10.6f %10.8f "
                "%12.6f %12.6f %12.6f %12.6f "
                "%10.8f " /* 10 columns */
                "%10.6f %10.8f %10.6f %10.6f "
                "%10.6f %10.8f %10.8f "
                "%12.6f %12.6f %12.6f %12.6f " /* 21 columns */
                "%10.6f %10.6f %10.6f %10.8f %10.6f %10.6f\n",
           (double) rho,
           (double) mu0, (double) mu0a, (double) mu0b, (double) mu0br,
           (double) ((mu0 - mu1)/delta),
           (double) ((mu0a - mu1a)/delta),
           (double) ((mu0b - mu1b)/delta),
           (double) dmu0,
           (double) s,
           (double) mu0c, (double) mu0d, (double) mu0e, (double) mueth,
           (double) mu0f, (double) mufr1, (double) mufr2,
           (double) ((mu0c - mu1c)/delta),
           (double) ((mu0d - mu1d)/delta),
           (double) ((mu0e - mu1e)/delta),
           (double) ((mu0f - mu1f)/delta),
           (double) pres0, (double) pres0a, (double) pres0b, (double) rpres0,
           (double) pres0c, (double) pres0d);
  }

  sphr_close(sphr);
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
  free(DCr);
  free(DTr);
  free(DCk);
  free(DTk);

  free(Br);
  free(dBr);
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
  return 0;
}

