/* fix corruptions in the snapshot files */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"
#include "ieutil.h"
#include "gsl/gsl_sf_bessel.h"



#ifdef D
int dim = D;
#else
int dim = 2;
#endif

int nmax = 10;
double rmax = 0;
xdouble Rmax = 0;
int numpt = 1024;
int dohnc = 0;
int singer = 0;
int ring = 0;
int mkcorr = 0;
int verbose = 0;
char *fnvir = NULL;
int snapshot = 0;



static void doargs(int argc, char **argv)
{
  argopt_t *ao = argopt_open(0);

  ao->desc = "fix corruptions in the snapshot files";
  argopt_add(ao, "-D", "%d", &dim, "dimension");
  argopt_add(ao, "-n", "%d", &nmax, "maximal order");
  argopt_add(ao, "-r", "%lf", &rmax, "rmax (flexible)");
  argopt_add(ao, "-R", "%" XDBLSCNF "f", &Rmax, "rmax (fixed)");
  argopt_add(ao, "-M", "%d", &numpt, "number of points along r");
  argopt_add(ao, "--hnc", "%b", &dohnc, "use the hypernetted chain approximation");
  argopt_add(ao, "--ring", "%b", &ring, "use the ring-sum formula");
  argopt_add(ao, "--sing", "%b", &singer, "use the Singer-Chandler formula for HNC");
  argopt_add(ao, "--corr", "%b", &mkcorr, "try to correct HNC");
  argopt_add(ao, "-o", NULL, &fnvir, "output virial coefficient");
  argopt_add(ao, "-s", "%b", &snapshot, "save intermediate snapshots");
  argopt_add(ao, "-v", "%b", &verbose, "be verbose");
  argopt_addhelp(ao, "--help");
  argopt_parse(ao, argc, argv);
  if ( rmax <= 0 ) rmax = nmax + 2;
  if ( mkcorr ) singer = ring = 0;
  if ( singer ) ring = 1;
  if ( verbose ) argopt_dump(ao);
  argopt_close(ao);
}



static int fixarr(int l, int l0, const char *fn,
   int npt, xdouble *arr)
{
  char *p, fnbak[80];
  int i;
  FILE *fp, *fq;

  if ( l == l0 ) return 0;
  die_if ( l > l0, "cannot fix %s, l %d > %d\n", fn, l, l0);

  /* backup the original file */
  strcpy(fnbak, fn);
  p = strrchr(fnbak, '.');
  if ( p == NULL ) p = fnbak + strlen(fn);
  sprintf(p, "l%d.dat", l0);
  i = rename(fn, fnbak);
  printf("moved %s to %s, return code %d\n", fn, fnbak, i);

  /* copy the first l arrays to the new file */
  xfopen(fp, fnbak, "rb", exit(-1));
  xfopen(fq, fn, "wb", exit(-1));
  for ( i = 0; i < l; i++ ) {
    SNAPSHOT_READARR(fp, fnbak, arr, npt);
    SNAPSHOT_WRITEARR(fq, fn, arr, npt);
  }
  printf("wrote new %s with %d arrays of size %d\n", fn, l, npt);
  fclose(fp);
  fclose(fq);
  return 0;
}



static void fixit(int nmax, int npt, xdouble rmax)
{
  xdouble *ckl = NULL, *tkl = NULL, *crl = NULL, *trl = NULL, *yrl = NULL;
  xdouble **ck = NULL, **tk = NULL, **cr = NULL, **tr = NULL, **yr = NULL;
  int lck, ltk, lcr, ltr, l = -1, l0;

  MAKE1DARR(ckl,  npt);
  MAKE1DARR(tkl,  npt);
  MAKE1DARR(crl,  npt);
  MAKE1DARR(trl,  npt);
  if ( dohnc || mkcorr ) {
    MAKE1DARR(yrl,  npt);
  }

  printf("npt = %d\n", npt);

  l0 = -1;
  SNAPSHOT_OPENF(lck, l0, fnck, ck, 0, ckl);
  l = lck;
  printf("loaded ck from %s, l %d\n", fnck, lck);

  l0 = -1;
  SNAPSHOT_OPENF(ltk, l0, fntk, tk, 0, tkl);
  die_if (ltk > lck, "ltk %d > lck %d\n", ltk, lck);
  if (ltk < lck) l = ltk;
  printf("loaded tk from %s, l %d\n", fntk, ltk);

  l0 = -1;
  SNAPSHOT_OPENF(lcr, l0, fncr, cr, 0, crl);
  die_if (lcr > lck, "lcr %d > lck %d\n", lcr, lck);
  if (lcr < lck) l = lcr;
  printf("loaded cr from %s, l %d\n", fncr, lcr);

  l0 = -1;
  SNAPSHOT_OPENF(ltr, l0, fntr, tr, 0, trl);
  die_if (ltr > lck, "ltr %d > lck %d\n", ltr, lck);
  if (ltr < lck) l = ltr;
  printf("loaded tr from %s, l %d\n", fntr, ltr);

  /* try to load yrl */
  SNAPSHOT_MKSTEM(fnyrstem, yr);
  if ( yrl != NULL ) {
    int n1 = 0, n2 = nmax - 1;
    while ( l > 0 ) {
      sprintf(fnyr, "%sl%d.dat", fnyrstem, l);
      if ( fexists(fnyr) ) break;
      l--;
    }
    die_if ( l == 0, "no yr for %s\n", fnyrstem);
    if ( snapshot_loadf(&n1, &n2, nmax, fnyr, npt, yr, 0, yrl) != 0 ) {
      fprintf(stderr, "cannot load yr(%d) of %d arrays from %s\n",
         l, n2, fnyr);
      exit(-1);
    }
    printf("loaded yr from %s, l %d\n", fnyr, l);
  }

  fixarr(l, lck, fnck, npt, ckl);
  fixarr(l, ltk, fntk, npt, tkl);
  fixarr(l, lcr, fncr, npt, crl);
  fixarr(l, ltr, fntr, npt, trl);

  FREE1DARR(ckl,  npt);
  FREE1DARR(tkl,  npt);
  FREE1DARR(crl,  npt);
  FREE1DARR(trl,  npt);
  if (yrl != NULL) FREE1DARR(yrl,  npt);
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
  printf("Rmax %g\n", (double) Rmax);
  fixit(nmax, numpt, Rmax);
  return 0;
}
