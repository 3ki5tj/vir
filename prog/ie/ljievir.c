/* Computing the virial coefficients of the Lennard-Jones fluid
 * by the PY or HNC integral equations
 * For the normal double precision
 *  gcc ljie3d.c -lfftw3
 * Or for the long double precision
 *  gcc -DLDBL ljie3d.c -lfftw3l
 * To disable FFTW
 *  gcc -DNOFFTW ljie3dfftw.c -lm
 * */
#include <stdio.h>
#include <math.h>
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "fftx.h"
#include "ieutil.h"



int dim = D;
int nmax = 10;
xdouble T = 1;
xdouble beta = 1;
xdouble rmax = 0;
int numpt = 4096;
int ffttype = 1;
int dohnc = 0;
int ring = 0;
int singer = 0;
int mkcorr = 0;
int fast = 0;
int verbose = 0;
char *fnvir = NULL;
char *fncrtr = NULL;

xdouble hncamp = 0, hncq = 1, hncalpha = -1;
xdouble shift = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);
  ao->desc = "computing the virial coefficients from the PY/HNC closure for the 3D hard-sphere fluid";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-T", "%" XDBLSCNF "f", &T, "temperature");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &rmax, "maximal r");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "-t", "%d", &ffttype, "FFT type");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted chain approximation");
  argopt_add(ao, "-q", "%" XDBLSCNF "f", &hncq,  "q of the hypernetted chain approximation, y(r) = a0 exp(q t(r))");
  argopt_add(ao, "-a", "%" XDBLSCNF "f", &hncamp, "a0 of the hypernetted chain approximation, 0: to be set as 1/q");
  argopt_add(ao, "-A", "%" XDBLSCNF "f", &hncalpha, "Roger-Young switch function 1 - exp(-alpha*r)");
  argopt_add(ao, "-c", "%" XDBLSCNF "f", &shift, "shift of t(r) in computing Bv");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--fast", "%b", &fast, "save tk and yr to accelerate the calculation");
  argopt_add(ao, "--corr", "%b", &mkcorr, "correct the closure");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "--crtr", NULL, &fncrtr, "file name of c(r) and t(r)");
  argopt_add(ao, "-v", "%+", &verbose, "be verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax <= 0 ) rmax = (nmax + 2)*2.0;
  beta = 1/T;
  if ( hncq > 0 && hncamp <= 0 )
    hncamp = 1/hncq; /* set amp automatically = 1/q */
  if ( mkcorr ) singer = ring = 0;
  else if ( FABS(hncamp - 1) > 1e-6 || FABS(hncq - 1) > 1e-6 || hncalpha >= 0 )
    dohnc = 1; /* turn on HNC, if necessary */
  if ( singer ) ring = 1;
  printf("rmax = %f, T %lf, HNC %d\n", (double) rmax, (double) T, dohnc);
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



/* return the potential phi(r), and -r*phi'(r)*/
static xdouble pot(xdouble r, xdouble *ndphir)
{
  xdouble invr6 = 1/(r*r);
  invr6 = invr6*invr6*invr6;
  *ndphir = invr6*(48*invr6 - 24);
  return 4*invr6*(invr6 - 1);
/*
  *ndphir = 2*r*r/(EXP(r*r) - 1);
  return -log(1-EXP(-r*r));
*/
/* for r^(-12)th potential */
/*
  *ndphir = 12*invr6*invr6;
  return invr6*invr6;
*/
}


/* save the header for the virial file */
__inline static char *LJsavevirhead(const char *fn,
    int dim, int l0, int nmax, xdouble T, int dohnc, int mkcorr, int npt,
    xdouble rmax, clock_t inittime)
{
  FILE *fp;
  static char fndef[256];

  if ( fn == NULL ) {
    sprintf(fndef, "LJBn%s%sD%dn%dT%gR%.0fM%d%s.dat",
        dohnc ? "HNC" : "PY", mkcorr ? "c" : "",
        dim, nmax, (double) T, (double) rmax, npt, STRPREC);
    fn = fndef;
  }
  xfopen(fp, fn, (l0 == 1) ? "w" : "a", return NULL);
  fprintf(fp, "# %s %s %d %.14f %d | n Bc Bv Bm [Bh Br | corr] | %.3fs\n",
      dohnc ? "HNC" : "PY", mkcorr ? "corr" : "",
      nmax, (double) rmax, npt, (double) inittime / CLOCKS_PER_SEC);
  fclose(fp);
  return (char *) fn;
}




/* compute the virial coefficients from the Percus-Yevick closure */
static int intgeq(int nmax, int npt, xdouble rmax, int ffttype, int dohnc)
{
  xdouble Bc, Bv, Bm = 0, Bh = 0, Br = 0, B2, B2tail = 0, fcorr;
  xdouble *fr, *rdfr, *swr = NULL, *fr1 = NULL, *rdfr1 = NULL, *phi2 = NULL;
  xdouble *crl, *trl, **ck, *tkl, **cr = NULL, **tr = NULL, **tk = NULL;
  xdouble *yrl, **yr = NULL, *vc = NULL;
  xdouble rm;
  int i, l;
  clock_t t0 = clock(), t1;
  sphr_t *sphr;

  sphr = sphr_open(dim, npt, rmax, 0, ffttype);
  rm = POW(2, (xdouble) 1 / 6);

  MAKE1DARR(fr, npt);
  MAKE1DARR(rdfr, npt);
  MAKE1DARR(fr1, npt);
  MAKE1DARR(phi2, npt);
  MAKE1DARR(rdfr1, npt);
  MAKE1DARR(swr, npt);
  MAKE1DARR(crl, npt);
  MAKE1DARR(trl, npt);
  MAKE1DARR(yrl, npt);
  MAKE1DARR(tkl, npt)
  MAKE2DARR(ck, nmax - 1, npt)
  if ( fast ) {
    MAKE2DARR(tk, nmax - 1, npt);
  }

  for ( i = 0; i < npt; i++ ) { /* compute f(r) = exp(-beta u(r)) - 1 */
    xdouble nrdphi, phi1, xp1;
    xdouble phi = pot(sphr->ri[i], &nrdphi);
    if ( sphr->ri[i] < rm ) {
      phi1 = phi + 1;
      phi2[i] = -1;
    } else {
      phi1 = 0;
      phi2[i] = phi;
      crl[i] = -beta * phi2[i];
    }
    xp1 = EXP(-beta * phi1);
    fr1[i] = xp1 - 1;
    rdfr1[i] = beta * nrdphi * xp1;
    fr[i] = EXP(-beta * phi) - 1;
    rdfr[i] = beta * nrdphi * (1 + fr[i]);
    if ( dohnc ) {
      crl[i] = xp1 * EXP(-beta * phi2[i]) - 1;
    } else {
      crl[i] = fr1[i] - (fr1[i] + 1) * beta * phi2[i];
    }
  }

  B2 = -integr(npt, fr, sphr->rDm1)/2;
  B2tail = -PI*8*pow_si(rmax, -3)/3;
  B2 += B2tail;

  printf("beta %g, dr %f, dm %d, rmax %f, ffttype %d, HNC %d, B2 %g %g\n",
      (double) beta, (double) sphr->dr, sphr->dm, (double) sphr->rmax,
      ffttype, dohnc, (double) B2, (double) B2tail);

  if ( singer ) {
    MAKE2DARR(cr, nmax - 1, npt);
    COPY1DARR(cr[0], fr, npt); /* TODO: fix this */
  }

  MAKE2DARR(tr, nmax - 1, npt);

  if ( dohnc || mkcorr ) {
    MAKE2DARR(yr, nmax - 1, npt);
    for ( i = 0; i < npt; i++ ) {
      swr[i] = 1;
      if ( hncalpha >= 0 )
        swr[i] = 1 - EXP( -hncalpha * sphr->ri[i] );
    }
    for ( i = 0; i < npt; i++ ) {
      yr[0][i] = hncamp / swr[i] * EXP( -beta * phi2[i] );
    }
    if ( mkcorr ) {
      MAKE1DARR(vc, npt);
    }
  }

  t1 = clock();
  fnvir = LJsavevirhead(fnvir, 3, 1, nmax, T,
      dohnc, mkcorr, npt, rmax, t1 - t0);

  for ( l = 1; l < nmax - 1; l++ ) {
    /* c(r) --> c(k) for the previous l */
    sphr_r2k(sphr, crl, ck[l-1]);

    /* compute the ring sum based on ck */
    if ( ring ) {
      Bh = get_BhBrk(l, npt, dohnc, ck, sphr->kDm1, sphr->ki, &Br);
    }

    /* compute t_l(k) from c_0(k), ... c_{l-1}(k) */
    if ( tk != NULL ) { /* fast version requires tk[] */
      get_tk_oz_fast(l, npt, ck, tk);
      COPY1DARR(tkl, tk[l], npt); /* tkl = tk[l] */
    } else { /* slow version requires ck[] only */
      get_tk_oz(l, npt, ck, tkl);
    }

    /* t_l(k) --> t_l(r) */
    sphr_k2r(sphr, tkl, trl);

    if ( tr != NULL ) {
      COPY1DARR(tr[l], trl, npt); /* tr[l] = trl; */
    }

    /* compute the cavity function y(r) */
    if ( yr != NULL ) {
      get_yr_hnc_fast(l, npt, yrl, yr, tr);
    } else {
      get_yr_hnc(l, npt, yrl, tr);
    }

    if ( mkcorr ) { /* construct the correction function */
      for ( i = 0; i < npt; i++ )
        vc[i] = yrl[i] - trl[i];
    }

    if ( dohnc ) {
      /* hypernetted-chain approximation: c(r) = (f(r) + 1) y(r) - (1 + t(r)) */
      for ( i = 0; i < npt; i++ )
        crl[i] = (fr1[i] + 1) * yrl[i] - trl[i];
      Bv = integr2(npt, rdfr1, yrl, sphr->rDm1) / (2*dim);
    } else {
      /* MSA approximation */
      for ( i = 0; i < npt; i++ )
        crl[i] = fr1[i] * trl[i];
      Bv = integr2(npt, rdfr1, trl, sphr->rDm1) / (2*dim);
    }

    /* B_{l+2}^c = -[1/(l+2)] Int c_l(r) 4 pi r^2 dr */
    Bc = -integr(npt, crl, sphr->rDm1) / (l + 2);

    if ( cr != NULL ) {
      COPY1DARR(cr[l], crl, npt); /* cr[l] = crl */
      if ( dohnc ) {
        Bm = get_Bm_singer(l, npt, cr, tr, sphr->rDm1);
        Bh += get_Bh_singer(l, npt, cr, tr, sphr->rDm1);
      } else {
        Bm = get_Bx_py(l, npt, cr, tr, sphr->rDm1);
      }
    } else {
      Bm = Bh = 0;
    }

    if ( mkcorr ) {
      if (l > 1) Bv *= 1 + shift;
      Bm = get_corr1(l, npt, crl, fr1, rdfr1, sphr->rDm1, 3,
          vc, 1, &Bc, &Bv, &fcorr);
    }

    savevir(fnvir, 3, l+2, Bc, Bv, Bm, Bh, Br, 0, B2, mkcorr, fcorr);
    savecrtr(fncrtr, l, npt, sphr->ri, crl, trl, vc, yrl);
  }
  savevirtail(fnvir, clock() - t1);

  sphr_close(sphr);
  FREE1DARR(fr, npt);
  FREE1DARR(rdfr, npt);
  FREE1DARR(fr1, npt);
  FREE1DARR(phi2, npt);
  FREE1DARR(rdfr1, npt);
  FREE1DARR(swr, npt);
  FREE1DARR(crl, npt);
  FREE1DARR(trl, npt);
  FREE1DARR(yrl, npt);
  FREE1DARR(tkl, npt);
  FREE2DARR(ck, nmax - 1, npt);
  FREE2DARR(tk, nmax - 1, npt);
  FREE2DARR(cr, nmax - 1, npt);
  FREE2DARR(tr, nmax - 1, npt);
  FREE2DARR(yr, nmax - 1, npt);
  FREE1DARR(vc, npt);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  intgeq(nmax, numpt, rmax, ffttype, dohnc);
  return 0;
}

