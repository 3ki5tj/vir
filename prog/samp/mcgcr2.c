/* grand-canonical simulation to compute the partition function
 * improved version for higher dimensions */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dgrjw.h"
#include "mcutil.h"



int nmin = 2; /* the minimal order of virial coefficients */
int nmax = DG_NMAX; /* the maximal order of virial coefficients */
int mtiers = 2;
real mcamp = 1.5f;
int nequil = 10000000;
double nsteps = 10000000;
double ratn = 0.5;
int gaussdisp = 0;
real rc0 = -1;
real sr0 = 1;
double mindata = 100; /* minimal # of data points to make update */
char *fninp = NULL;
char *fnZrr = NULL;
char *fnZr = NULL;

int nstcom = 10000;
int nsted = 10;
int nstcs = 100;
int nstfb = -1;
/* maximal number of extra edges to invoke fb calculation
 * reference T60 timing data (no clique separator):
 * nedxmax  time per evaluation (ms)
 *  10         15
 *  11         50
 *  12        150
 * */
int nedxmax = -1;

int updzr = 0; /* update rc */
int neql = 0; /* rounds of equilibration (in which Z is updated) */
int nstsave = 1000000000; /* interval of saving data */

int restart = 0; /* restartable simulation */
int bsim0 = 0; /* first simulation */



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int m1 = 1;

  ao = argopt_open(0);
  argopt_add(ao, "-z",    "%d",   &nmin,    "minimal order of virial coefficients");
  argopt_add(ao, "-n",    "%d",   &nmax,    "maximal order of virial coefficients");
  argopt_add(ao, "-m",    "%d",   &m1,      "number of intermediate ensembles per order of virial coefficients");
  argopt_add(ao, "-1",    "%lf",  &nsteps,  "number of simulation steps");
  argopt_add(ao, "-a",    "%r",   &mcamp,   "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-r",    "%lf",  &ratn,    "rate of particle moves");
  argopt_add(ao, "-c",    "%r",   &rc0,     "radius of particle insersion");
  argopt_add(ao, "-s",    "%r",   &sr0,     "default scaling of distance between the constrained vertex pair");
  argopt_add(ao, "-i",    NULL,   &fninp,   "input file");
  argopt_add(ao, "-o",    NULL,   &fnZrr,   "output file");
  argopt_add(ao, "-O",    NULL,   &fnZr,    "output file (compact)");
  argopt_add(ao, "-G",    "%d",   &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-A",    "%d",   &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-F",    "%d",   &nstfb,   "interval of computing fb");
  argopt_add(ao, "-X",    "%d",   &nedxmax, "maximal number of extra edges for computing fb");
  argopt_add(ao, "-Q",    "%lf",  &mindata, "minimal number of data points to warrant parameter update");
  argopt_add(ao, "-q",    "%d",   &nstsave, "interval of saving data");
  argopt_add(ao, "-E",    "%d",   &neql,    "number of equilibration rounds");
  argopt_add(ao, "-0",    "%d",   &nequil,  "number of equilibration steps per round");
  argopt_add(ao, "-M",    "%d",   &nstcom,  "interval of centering the structure");
  argopt_add(ao, "-R",    "%b",   &restart, "run a restartable simulation");
  argopt_add(ao, "-B",    "%b",   &bsim0,   "start a restartable simulation");
  argopt_add(ao, "-C",    "%b",   &updzr,   "update Zr for pure states");

  argopt_parse(ao, argc, argv);

  die_if (nmax > DG_NMAX, "too many %d atoms >= %d\n", nmax, DG_NMAX);

  mtiers = m1 + 1;
  die_if (mtiers <= 1, "mtiers %d should be at least 2\n", mtiers);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  if (rc0 <= 0) {
    rc0 = (real) (4. / sqrt(D));
    if (rc0 > 1) rc0 = 1;
  }

  /* in higher dimensions compute fb more often */
  if (D > 20) nstcs = nsted;
  if (nstfb < 0) nstfb = nstcs;

  if (nedxmax < 0) {
    int nedarr[] = {0, 9, 9, 9, 9, 9, 9, 11, 11, 12, 12, 12, 13, 13};
    nedxmax = (D < 14) ? nedarr[D] : 14;
  }

  argopt_dump(ao);
  argopt_close(ao);
}



/* object-oriented definition and functions */



#define GCX_PURE 0



typedef struct {
  int nmin, nmax;
  int m; /* number of ensembles for each n */
  int ens0; /* starting ensemble index */
  int nens; /* number of ensembles */
  int *type; /* GCX_PURE == 0 for pure states */
  int *n; /* number of vertices */
  double *Z; /* partition function */
  double *hist; /* histogram */
  double (*ndown)[2]; /* number of upward (vertex adding) moves */
  double (*nup)[2]; /* number of downward (vertex removing) moves */
  double *Zr; /* Z[iens] / Z[iens - 1] */
  real *rc; /* volume of the insertion */
  real *sr; /* scaling between the special pair of vertices */
  double *Zr1; /* alternative buffer for Zr */
  real *rc1, *sr1; /* alternative buffer for rc and sr */
  double nedg[DG_NMAX + 1][2]; /* number of edges */
  double ncsp[DG_NMAX + 1][2]; /* no clique separator */
  double fbsm[DG_NMAX + 1][3]; /* hard-sphere weight */
  double B[DG_NMAX + 1]; /* virial coefficients */
} gc_t;



/* clear data */
static void gc_cleardata(gc_t *gc)
{
  int i;

  for (i = 0; i < gc->nens; i++) {
    gc->hist[i] = 0;
    gc->nup[i][0] = gc->ndown[i][0] = 1e-20;
    gc->nup[i][1] = gc->ndown[i][1] = 0;
  }
  for (i = 0; i <= gc->nmax; i++) {
    gc->nedg[i][0] = 1e-20;
    gc->nedg[i][1] = 0;
    gc->ncsp[i][0] = 1e-20;
    gc->ncsp[i][1] = 0;
    gc->fbsm[i][0] = 1e-20;
    gc->fbsm[i][1] = gc->fbsm[i][2] = 0;
  }
}



INLINE gc_t *gc_open(int nmin, int nmax, int m,
    real rc0, real sr0)
{
  gc_t *gc;
  int i, n, nensmax;
  double rc, sr;

  xnew(gc, 1);
  gc->nmin = nmin;
  gc->nmax = nmax;
  gc->m = m;
  gc->ens0 = nmin * m;
  gc->nens = 1 + nmax * m;
  nensmax = 1 + DG_NMAX * m;
  xnew(gc->type, nensmax);
  xnew(gc->n, nensmax);
  xnew(gc->Z, nensmax);
  xnew(gc->hist, nensmax);
  xnew(gc->nup, nensmax);
  xnew(gc->ndown, nensmax);
  xnew(gc->Zr, nensmax);
  xnew(gc->rc, nensmax);
  xnew(gc->sr, nensmax);
  xnew(gc->Zr1, nensmax);
  xnew(gc->rc1, nensmax);
  xnew(gc->sr1, nensmax);
  /* temporary setting */
  if (m > 2 && rc0 >= 1) rc0 = 0.8;
  for (i = 0; i < nensmax; i++) {
    gc->n[i] = (i + m - 1) / m;
    /* type == 0 mean GCX_PURE */
    gc->type[i] = (i % m);
  }

  for (i = 0; i < nensmax; i++) {
    n = gc->n[i];
    /* guess rc */
    if ( gc->type[i] == GCX_PURE ) {
      rc = 1;
    } else {
      /* rc = pow(rc0, 1.*(m - i % m)/(m - 1)); */
      /* linear interpolation, including the boundary value 1 */
      rc = rc0 + (1. - rc0) * (i % m - 1) / (m - 1);
    }
    gc->rc[i] = (real) rc;

    if (gc->type[i] == GCX_PURE)
      gc->Zr[i] = 1.0/n; /* inverse number of nonarticulated pairs */
    else gc->Zr[i] = 1;
    gc->Z[i] = 1; /* absolute partition function */
  }

  /* sr requires a separate loop for it needs rc of a higher ensemble */
  for (i = 0; i < nensmax; i++) {
    if ( i == nensmax - 1 || gc->type[i] == GCX_PURE ) {
      sr = 1;
    } else if ( sr0 > 1 && i > 0 && gc->type[i - 1] == GCX_PURE ) {
      /* the given value is only used for the first rc */
      sr = sr0;
    } else { /* guess from the ratios of the successive radii */
      if (gc->rc[i] > 1) sr = 1;
      else sr = gc->rc[i + 1] / gc->rc[i];
    }
    gc->sr[i] = (real) sr;
  }

  for (i = 0; i <= DG_NMAX; i++) gc->B[i] = 1;
  gc_cleardata(gc);
  return gc;
}



INLINE void gc_close(gc_t *gc)
{
  free(gc->type);
  free(gc->n);
  free(gc->Z);
  free(gc->hist);
  free(gc->nup);
  free(gc->ndown);
  free(gc->Zr);
  free(gc->rc);
  free(gc->sr);
  free(gc->Zr1);
  free(gc->rc1);
  free(gc->sr1);
  free(gc);
}



/* accumulate data */
static void gc_accumdata(gc_t *gc, const dg_t *g, double t,
    int nstcs, int nstfb)
{
  int ned = -1, n = g->n, ncs, sc, err, fb = 0;
  int nlookup = DGMAP_NMAX;
  static int degs[DG_NMAX];

  /* we will always compute fb if a lookup table is available
   * i.e., for a small n, but the n = 8 case costs 30s-70s,
   * which is a bit long for a short simulation */
  if (nsteps < 5e8) nlookup = 7;

  ned = dg_degs(g, degs);
  gc->nedg[n][0] += 1;
  gc->nedg[n][1] += ned;

  if (n <= nlookup || rnd0() * nstcs / nsted < 1) {
    /* check if the graph has a clique separator */
    if (n <= 3) {
      int fbarr[4] = {1, 1, -1, -1};
      err = 0;
      ncs = 1;
      fb = fbarr[n];
    } else if (n <= nlookup) { /* lookup table */
      code_t code;
      unqid_t uid;
      dgmap_t *m = dgmap_ + n;

      dg_encode(g, &code);
      ncs = (dg_ncsep_lookuplow(g, code) == 0);
      dgmap_init(m, n);
      uid = m->map[code];
      fb = dg_hsfb_lookuplow(n, uid);
      err = 0;
    } else {
      /* this function implicitly computes the clique separator
       * with very small overhead */
      sc = dg_rhsc_spec0(g, 0, &ned, degs, &err);
      if (err == 0) {
        ncs = (sc != 0);
        fb = DG_SC2FB(sc, ned);
      } else { /* no clique separator */
        ncs = 1;
      }
    }
    gc->ncsp[n][0] += 1;
    gc->ncsp[n][1] += ncs;

    if (n <= nlookup || rnd0() * nstfb/nstcs < 1) {
      /* compute fb, if it is cheap */
      if ( err ) { /* if dg_rhsc_spec0() fails, no clique separator */
        if (ned > n + nedxmax && n > nedxmax + 2) {
          fprintf(stderr, "t %.6e/%g D %d assume fb = 0 n %d ned %d\n",
              t, nsteps, D, n, ned);
          fb = 0;
        } else {
          fb = dg_hsfb_mixed0(g, 1, &ned, degs);
        }
      }
      gc->fbsm[n][0] += 1;
      gc->fbsm[n][1] += fb;
      gc->fbsm[n][2] += abs(fb);
    }
  }
}



/* update parameters */
static void gc_update(gc_t *gc, double mindata, int updzr,
    double *Zr, real *rc, real *sr)
{
  int i;
  double r;

  if (Zr != gc->Zr) memcpy(Zr, gc->Zr, sizeof(Zr[0]) * gc->nens);
  if (rc != gc->rc) memcpy(rc, gc->rc, sizeof(rc[0]) * gc->nens);
  if (sr != gc->sr) memcpy(sr, gc->sr, sizeof(sr[0]) * gc->nens);
  for (i = gc->ens0 + 1; i < gc->nens; i++) {
    /* make sure enough data points and not to change the exact data */
    if (gc->ndown[i][1] >= mindata && gc->nup[i-1][1] >= mindata) {
      r = (gc->nup[i-1][1] / gc->nup[i-1][0])
        / (gc->ndown[i][1] / gc->ndown[i][0]);
      if (gc->type[i - 1] == GCX_PURE) {
        Zr[i] *= r;
      } else {
        sr[i - 1] *= pow(r, 1./D);
      }
    }
  }
}



/* compute the partition function and virial coefficients */
static void gc_computeZ(gc_t *gc,
    const double *Zr, const real *rc, const real *sr)
{
  double fac = 1, z = 1, fbav;
  int n, i;

  if (Zr == NULL) Zr = gc->Zr;
  if (rc == NULL) rc = gc->rc;
  if (sr == NULL) sr = gc->sr;
  for (i = 1; i < gc->nens; i++) {
    n = gc->n[i];
    /* compute the partition function */
    z *= Zr[i];
    if (gc->type[i - 1] != GCX_PURE) {
      z *= pow(sr[i - 1], D);
    } else {
      z *= pow(rc[i], D);
    }
    //printf("compute Z, i %d, n %d, rc %.7f, sr %.7f, Z %.6e\n", i, n, rc[i], sr[i - 1], x);
    /* use the the exact result if available */
    if (gc->type[i] == GCX_PURE) {
      z *= n * (n - 1) * .5; /* pair generation bias */
      if (n == 3 && nmin < 3)
        printf("Z3 %g vs %g(ref)\n", z, Z3rat(D));
      if (n == 4)
        printf("Z4 %g vs %g(ref, %s)\n", z, Z4rat(D), (D > 12) ? "approx." : "exact");
      if (n <= 2) z = 1;
      else if (n == 3) z = Z3rat(D);
      else if (n == 4 && D <= 12) z = Z4rat(D);
    }
    gc->Z[i] = z;
    if (gc->type[i] == GCX_PURE) {
      if (n >= 2) fac *= 2./n;
      fbav = gc->fbsm[n][1] / gc->fbsm[n][0];
      /* Bn/B2^(n-1) = (1 - n) 2^(n - 1) / n! Zn/Z2^(n-1) < fb > */
      gc->B[n] = (1. - n) * fac * z * fbav;
      //printf("B[%d] = %g\n", n, gc->B[n]);
    }
  }
}



/* print out a summary */
static void gc_print(const gc_t *gc, int compact,
    const double *Zr, const real *rc, const real *sr)
{
  int i;

  if (Zr == NULL) Zr = gc->Zr;
  if (rc == NULL) rc = gc->rc;
  if (sr == NULL) sr = gc->sr;
  for (i = gc->ens0; i < gc->nens; i++) {
    printf("%3d %2d %9.6f %8.6f %8.6f ",
       i, gc->n[i], Zr[i], rc[i], sr[i]);
    if ( !compact ) {
      printf("%12.0f %.5f %.5f ", gc->hist[i],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0]);
      if (gc->type[i] == GCX_PURE) {
        int n = gc->n[i];
        printf("%6.3f %8.6f %+9.7f %+.5e",
          gc->nedg[n][1] / gc->nedg[n][0],
          gc->ncsp[n][1] / gc->ncsp[n][0],
          gc->fbsm[n][1] / gc->fbsm[n][0],
          gc->B[n]);
      }
    }
    printf("\n");
  }
}



#define mkfnZrdef(fn, fndef, d, m, n) \
  if (fn == NULL) { sprintf(fndef, "ZrD%dr%dn%d.dat", d, m, n); fn = fndef; }



/* save a compact (pure-state) Zr file, mimic the output of mcgc.c */
static int gc_saveZr(gc_t *gc, const char *fn)
{
  FILE *fp;
  int i, n;
  char fndef[64];
  double Zr = 1, prevZ = 1, tot = 0;

  mkfnZrdef(fn, fndef, D, gc->m - 1, gc->nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d V2 %d\n", D, gc->nmax, nedxmax);
  for (i = 1; i < gc->nens; i++) {
    if (gc->type[i] != GCX_PURE) continue;
    n = gc->n[i];
    /* start with the third virial coefficient */
    if (n <= 2) continue;
    Zr = gc->Z[i] / prevZ;
    fprintf(fp, "%3d %18.14f %20.14e %14.0f %14.0f %14.0f %.14f %.14f "
        "%14.0f %14.0f %14.0f %18.14f %16.14f %16.14f %+17.14f %+20.14e\n",
        n, Zr, gc->Z[i], gc->hist[i],
        gc->nup[i-1][0], gc->ndown[i][0],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0],
        gc->nedg[n][0], gc->ncsp[n][0], gc->fbsm[n][0],
        gc->nedg[n][1] / gc->nedg[n][0], gc->ncsp[n][1] / gc->ncsp[n][0],
        gc->fbsm[n][2] / gc->fbsm[n][0], gc->fbsm[n][1] / gc->fbsm[n][0],
        gc->B[n]);
    tot += gc->hist[i];
    prevZ = gc->Z[i];
  }
  fclose(fp);
  printf("saved compact Zr to %s, total %g\n", fn, tot);
  return 0;
}



#define mkfnZrrdef(fn, fndef, d, m, n) \
  if (fn == NULL) { sprintf(fndef, "ZrrD%dr%dn%d.dat", d, m, n); fn = fndef; }



/* save all data to file */
static int gc_saveZrr(gc_t *gc, const char *fn,
    const double *Zr, const double *rc, const double *sr)
{
  FILE *fp;
  int i, n;
  char fndef[64];
  double tot = 0;

  mkfnZrrdef(fn, fndef, D, gc->m - 1, gc->nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#R %d %d %d %d %d %d V2 %d\n",
      D, gc->nens, gc->ens0, gc->nmin, gc->nmax, gc->m, nedxmax);
  for (i = 1; i < gc->nens; i++) {
    fprintf(fp, "%4d %d %3d %20.14e %20.14e %18.14f %18.14f "
        "%14.0f %14.0f %14.0f %.14f %.14f ",
        i, gc->type[i], gc->n[i], Zr[i], gc->Z[i], rc[i], sr[i],
        gc->hist[i], gc->nup[i-1][0], gc->ndown[i][0],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0]);
    if (gc->type[i] == GCX_PURE) {
      n = gc->n[i];
      fprintf(fp, "%14.0f %14.0f %14.0f ",
        gc->nedg[n][0], gc->ncsp[n][0], gc->fbsm[n][0]);
      fprintf(fp, "%18.14f %16.14f %16.14f %+17.14f %+20.14e",
        gc->nedg[n][1] / gc->nedg[n][0], gc->ncsp[n][1] / gc->ncsp[n][0],
        gc->fbsm[n][2] / gc->fbsm[n][0], gc->fbsm[n][1] / gc->fbsm[n][0],
        gc->B[n]);
    }
    fprintf(fp, "\n");
    tot += gc->hist[i];
  }
  fclose(fp);
  printf("saved Zrr to %s, tot %g\n", fn, tot);
  return 0;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int gc_load(gc_t *gc, const char *fn, int loadall)
{
  FILE *fp;
  int i = -1, nens, d = 0, n, tp, ens0 = 0, n0, m = 0, ver = 0, next;
  char s[512], fndef[64], sver[16], *p;

  mkfnZrrdef(fn, fndef, D, gc->m - 1, gc->nmax);
  xfopen(fp, fn, "r", return -1);

  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if ( s[0] != '#' || s[1] != 'R'
    || 7 != sscanf(s + 2, "%d%d%d%d%d%d%8s",
                   &d, &nens, &ens0, &n0, &n, &m, sver)
    || d != D || m != gc->m || ens0 != gc->ens0 ) {
    fprintf(stderr, "%s: D %d vs. %d, m %d vs %d, ens0 %d vs %d\n%s",
        fn, d, D, m, gc->m, ens0, gc->ens0, s);
    fclose(fp);
    return -1;
  }
  ver = atoi(sver[0] == 'V' ? sver + 1 : sver); /* strip V */
  if (loadall && ver <= 0) {
    fprintf(stderr, "%s version %s does not have restartable data\n", fn, sver);
    loadall = 0;
  }

  for (i = 1; i < nens; i++) {
    double Zr, Z, rc, sr;

    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if (7 != sscanf(s, "%d%d%d%lf%lf%lf%lf%n",
                       &d, &tp, &n, &Zr, &Z, &rc, &sr, &next)
        || i != d || tp != gc->type[i] || n != gc->n[i]) {
      fprintf(stderr, "%s ends on line %d\n%s", fn, i, s);
      break;
    }
    if (Zr >= 0) gc->Zr[i] = Zr;
    if (rc >= 0) gc->rc[i] = rc;
    if (sr >= 0) {
      if (ver >= 2) gc->sr[i] = sr;
      else gc->sr[i - 1] = sr; /* old version */
    }

    if (loadall) { /* try to get additional data */
      p = s + next;
      if (5 != sscanf(p, "%lf%lf%lf%lf%lf%n",
            &gc->hist[i], &gc->nup[i-1][0], &gc->ndown[i][0],
            &gc->nup[i-1][1], &gc->ndown[i][1], &next) ) {
        fprintf(stderr, "%s corrupted on line %d\n%s", fn, i, s);
        break;
      }
      if (gc->nup[i-1][0] <= 0) gc->nup[i-1][0] = 1e-20;
      if (gc->ndown[i][0] <= 0) gc->ndown[i][0] = 1e-20;
      gc->nup[i-1][1] *= gc->nup[i-1][0];
      gc->ndown[i][1] *= gc->ndown[i][0];
      if (gc->type[i] == GCX_PURE) {
        p += next;
        if (7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf",
              &gc->nedg[n][0], &gc->ncsp[n][0], &gc->fbsm[n][0],
              &gc->nedg[n][1], &gc->ncsp[n][1],
              &gc->fbsm[n][2], &gc->fbsm[n][1]) ) {
          fprintf(stderr, "%s corrupted on line %d\n%s", fn, i, s);
          break;
        }
        if (gc->nedg[n][0] <= 0) gc->nedg[n][0] = 1e-20;
        if (gc->ncsp[n][0] <= 0) gc->ncsp[n][0] = 1e-20;
        if (gc->fbsm[n][0] <= 0) gc->fbsm[n][0] = 1e-20;
        gc->nedg[n][1] *= gc->nedg[n][0];
        gc->ncsp[n][1] *= gc->ncsp[n][0];
        gc->fbsm[n][1] *= gc->fbsm[n][0];
        gc->fbsm[n][2] *= gc->fbsm[n][0];
      }
    }
  }
  fclose(fp);
  if (i >= nens) /* recompute the partition function */
    gc_computeZ(gc, gc->Zr, gc->rc, gc->sr);
  printf("loaded Zr from %s, version %d\n", fn, ver);
  return i;
}



/* single-vertex move with a distance restraint */
INLINE int bcrstep(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xi, real r2ij[][DG_NMAX], real r2i[],
    real rc, real amp, int gauss)
{
  int i, j, n = g->n;

  DISPRNDI(i, n, x, xi, amp, gauss);
  /* check if the distance constraint is satisfied */
  if (i == i0) { /* check xi and x[j0] */
    if ( rvn_dist2(xi, x[j0]) >= rc * rc )
      return 0;
  } else if (i == j0) { /* check xi and x[i0] */
    if ( rvn_dist2(xi, x[i0]) >= rc * rc )
      return 0;
  }
  UPDGRAPHR2(i, n, g, ng, x, xi, 1, r2i);
  /* the new graph needs to be biconnected and without
   * the articulation pair (i0, j0) */
  if ( dg_biconnected(ng)
    && dg_connectedvs(ng, mkbitsmask(n) ^ MKBIT(i0) ^ MKBIT(j0)) ) {
    dg_copy(g, ng);
    rvn_copy(x[i], xi);
    UPDR2(r2ij, r2i, n, i, j);
    return 1;
  }
  return 0;
}



/* attach a vertex to a random vertex i0, become a restrained state
 *  pure:
 *    bc(r^n) dr^n
 *  restrained:
 *    bc(R^{n+1}) step(Rc - R_{i0, n}) c(R^{n+1}\{i0, n}) dR^{n+1}
 *  transition state:
 *    bc(r^n) step(Rc - R_{i0, n}) bc(R^{n+1}) dR^{n+1}
 *  since bc(r^n) implies c(r^n\{i0}) == c(R^{n+1}\{i0, n})
 * and the transition probability is < bc(R^{n+1}) > */
INLINE int nmove_pureup2restrained(int *i0, dg_t *g,
    rvn_t *x, real r2ij[][DG_NMAX], real rc, double Zr)
{
  int i, n = g->n, deg = 0;

  /* the Zr step is cheaper, so we do it first */
  if (Zr < 1 && rnd0() >= Zr) return 0;

  *i0 = (int) (rnd0() * n); /* randomly choose a vertex as the root */
  /* xn = x0 + rc * rndball */
  rvn_inc( rvn_rndball(x[n], rc), x[*i0] );
  /* biconnectivity means to be linked to two vertices */
  for (i = 0; i < n; i++)
    if ( (r2ij[n][i] = r2ij[i][n] = rvn_dist2(x[i], x[n])) < 1 )
      deg++;
  if ( deg < 2 ) return 0;
  /* update the graph */
  g->n = n + 1;
  for (i = 0; i < n; i++) {
    if ( r2ij[n][i] < 1 ) dg_link(g, n, i);
    else dg_unlink(g, n, i);
  }
  return 1;
}



/* remove the one vertex out of the restrained bond (i0, j0),
 * and move from a mixed state down to a pure state
 *  restrained:
 *    bc(R^n) step(Rc - R_{i0, j0}) c(R^{n}\{i0,j0}) dR^n
 *  pure:
 *    bc(r^{n-1}) dr^{n-1}
 *  transition state:
 *    bc(r^{n-1}) step(Rc - R_{i0, j0}) bc(R^n) dR^n
 * and the transition probability is < bc(r^{n-1}) > */
INLINE int nmove_restraineddown2pure(int i0, int j0,
    dg_t *g, rvn_t *x, real r2ij[][DG_NMAX], double Zr)
{
  int i, n = g->n, j, k;

  if (Zr < 1 && rnd0() >= Zr) return 0;
  i = (rnd0() > 0.5) ? i0 : j0;
  if ( !dg_biconnectedvs(g, mkbitsmask(n) ^ MKBIT(i)) )
    return 0;
  dg_remove1(g, g, i);
  /* update x */
  for (j = i; j < n - 1; j++)
    rvn_copy(x[j], x[j + 1]);
  /* update the r2ij matrix */
  for (j = i; j < n - 1; j++) {
    for (k = 0; k < i; k++)
      r2ij[j][k] = r2ij[k][j] = r2ij[j + 1][k];
    for (k = i; k < j; k++)
      r2ij[j][k] = r2ij[k][j] = r2ij[j + 1][k + 1];
    r2ij[j][j] = r2ij[j + 1][j + 1];
  }
  return 1;
}



#define changenapair(i, j, g, r2ij, rc) \
    changenapair_metro(i, j, g, r2ij, rc)

/* change the distance-restrained pair */
INLINE int changenapair_simple(int *i, int *j, const dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int tmp;

  (void) g; (void) r2ij; (void) rc;
  tmp = *i, *i = *j, *j = tmp;
  return 1;
}



/* change the distance-restrained pair */
INLINE int changenapair_metro(int *i, int *j, const dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int ni, nj;

  ni = randpair(g->n, &nj);
  if ( !(ni == *i && nj == *j) && !(ni == *j && nj == *i)
    && dg_connectedvs(g, mkbitsmask(g->n) ^ MKBIT(ni) ^ MKBIT(nj))
    && r2ij[ni][nj] < rc * rc ) {
    *i = ni;
    *j = nj;
    return 1;
  }
  ni = *i, *i = *j, *j = ni;
  return 0;
}



/* scale the distance between x[i0] and x[j0] and test the biconnectivity */
INLINE int nmove_scale(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xj0, real r2ij[][DG_NMAX], real *r2i,
    real sr, real Zr, int extra)
{
  int k, n = g->n;

  if (Zr < 1 && rnd0() >= Zr) return 0;
  /* randomly choose i0 or j0 as the root */
  if (rnd0() < 0.5) { k = i0, i0 = j0, j0 = k; }
  /* xj0 = x[i0] + (x[j0] - x[i0]) * s */
  rvn_diff(xj0, x[j0], x[i0]);
  rvn_inc(rvn_smul(xj0, sr), x[i0]);
  /* make a new graph with the new xj0 */
  ng->n = g->n;
  dg_copy(ng, g);
  for (k = 0; k < n; k++) {
    /* update connection from j0 to others */
    if ( k == j0 ) continue;
    if ( (r2i[k] = rvn_dist2(x[k], xj0)) < 1 )
      dg_link(ng, k, j0);
    else
      dg_unlink(ng, k, j0);
  }
  if ( !dg_biconnected(ng) ) return 0;
  if ( extra && rnd0() >= 1.*Zr/dg_nedges(ng) ) return 0;
  dg_copy(g, ng);
  rvn_copy(x[j0], xj0);
  UPDR2(r2ij, r2i, n, j0, k);
  return 1;
}



#ifdef CHECK
/* check if everything is okay */
INLINE int check(const dg_t *g, rvn_t *x, real (*r2ij)[DG_NMAX],
    const gc_t *gc, int iens, int i0, int j0,
    double t, const char *msg)
{
  int i, j, n = g->n, err = 0, type = gc->type[iens];
  real r2, rc2 = gc->rc[iens] * gc->rc[iens];

  /* check the pair distances */
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      r2 = rvn_dist2(x[i], x[j]);
      if (r2ij[i][j] != r2ij[j][i]) {
        fprintf(stderr, "mat %d and %d, r2 %g vs %g\n",
            i, j, r2ij[i][j], r2ij[j][i]);
        err++;
      }
      if (fabs(r2ij[i][j] - r2) > 1e-6) {
        fprintf(stderr, "r2 %g (actual) vs %g r2ij\n",
            r2, r2ij[i][j]);
        err++;
      }
      if ( (r2 < 1 && !dg_linked(g, i, j))
        || (r2 >= 1 && dg_linked(g, i, j)) ) {
        fprintf(stderr, "bad linkage between %d, %d, r2 %g\n", i, j, r2);
        err++;
      }
    }
  }

  if ( !dg_biconnected(g) ) {
    fprintf(stderr, "diagram not biconnected\n");
    dg_print(g);
    err++;
  }

  if (type != 0) {
    if (i0 == j0 || i0 < 0 || j0 < 0 || i0 >= n || j0 >= n) {
      fprintf(stderr, "i0 %d vs j0 %d, n %d\n", i0, j0, n);
      err++;
    }
    if (r2ij[i0][j0] >= rc2) {
      fprintf(stderr, "broken constraint: i0 %d, j0 %d, n %d, "
          "r2 %g >= %g, rc = %g\n",
          i0, j0, n, r2ij[i0][j0], rc2, gc->rc[iens]);
      err++;
    }
    if ( !dg_connectedvs(g, mkbitsmask(g->n) ^ MKBIT(i0) ^ MKBIT(j0)) ) {
      fprintf(stderr, "the pair (%d, %d) is articulated\n", i0, j0);
      err++;
    }
  }
  die_if (err, "t %g, iens %d, type %d, n %d, message %s, err %d\n",
      t, iens, gc->type[iens], gc->n[i], msg, err);
  return err;
}
#endif



static void mcgcr(int nmin, int nmax, int mtiers, double nsteps,
    real mcamp, int neql, int nequil, int nstsave)
{
  double t, ctot = 0, cacc = 0, prtot = 0, pracc = 0;
  int iens, ensmax, ensmin, acc, it, ieql;
  int pi = 0, pj = 1; /* indices of the distance-restrained pair */
  gc_t *gc;
  dg_t *g, *ng, *g1, *ng1;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  real r2ij[DG_NMAX][DG_NMAX] = {{0}}, r2i[DG_NMAX] = {0};

  gc = gc_open(nmin, nmax, mtiers, rc0, sr0);
  gc_load(gc, fninp, restart && !bsim0);
  gc_print(gc, 1, NULL, NULL, NULL);
  ensmin = nmin * mtiers;
  ensmax = nmax * mtiers;
  iens = ensmin; /* lowest ensemble */

  g = dg_open(nmax);
  ng = dg_open(nmax);
  g1 = dg_open(nmax);
  ng1 = dg_open(nmax);
  /* initial n is nmin */
  calcr2ij(r2ij, x, nmin);
  mkgraphr2ij(g, r2ij, 1, nmin);
  printf("iens %d, ensmin %d, ensmax %d, n %d\n", iens, ensmin, ensmax, g->n);

  ieql = (neql > 0); /* rounds of equilibrations */

#ifdef CHECK
  fprintf(stderr, "checking code enabled\n");
#endif

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    //printf("t %g, iens %d, ensmin %d, ensmax %d, n %d\n", t, iens, ensmin, ensmax, g->n);
    die_if (g->n < nmin || g->n > nmax, "bad n %d, t %g, iens %d, mtiers %d, ensmax %d\n", g->n, t, iens, mtiers, ensmax);
    if (rnd0() < ratn) /* n-move, switching ensemble */
    {
      if (rnd0() < 0.5) { /* increase the ensemble index */
        if (iens >= ensmax) goto STEP_END;

        acc = 0;
        if (gc->type[iens] == GCX_PURE) { /* pure up to restrained */
          /*  bc(r^n) dr^n
           *    -->
           *  bc(R^{n+1}) step(Rc - R_{i0, n}) c(R^{n+1}\{i0, n}) dR^{n+1} */
          pj = g->n;
          acc = nmove_pureup2restrained(&pi, g, x, r2ij,
              gc->rc[iens + 1], 1. / gc->Zr[iens + 1]);
          //printf("t %g, n+ move 1, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gc, iens + acc, pi, pj, t, "nmove_pureup2restrained");
#endif
          die_if (g->n < nmin || g->n > nmax, "p2m BAD n %d, t %g, iens %d\n", g->n, t, iens);
        } else if (gc->type[iens + 1] == GCX_PURE) { /* restrained up to pure */
          /*  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           *    -->
           *  bc(R^n) dR^n */
          if ( r2ij[pi][pj] * dblsqr(gc->sr[iens]) < dblsqr(gc->rc[iens + 1]) )
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
              gc->sr[iens], 1. / gc->Zr[iens + 1], 1);
        } else { /* restrained up to a less restrained state */
          /*  c(r^n\{i0, j0}) step(rc - r_{i0, j0}) bc(r^n) dr^n
           *    -->
           *  c(R^n\{i0, j0}) step(Rc - R_{i0, j0}) bc(R^n) dR^n
           * where rc < Rc */
          if ( r2ij[pi][pj] * dblsqr(gc->sr[iens]) < dblsqr(gc->rc[iens + 1]) )
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
              gc->sr[iens], 1., 0);
          //printf("t %g, n+ move 2, acc %d, iens %d, n %d, n1 %d\n", t, acc, iens, g->n, g1->n);
#ifdef CHECK
          check(g, x, r2ij, gc, iens + acc, pi, pj, t, "nmove_scaleup");
#endif
        }
        gc->nup[iens][0] += 1;
        gc->nup[iens][1] += acc;
        if ( acc ) iens++;
        die_if (g->n < nmin || g->n > nmax, "BAD n %d, t %g, iens %d\n", g->n, t, iens);
      }
      else /* decrease the ensemble index */
      {
        if (iens <= ensmin) goto STEP_END;

        acc = 0;
        if (gc->type[iens - 1] == GCX_PURE) { /* restrained down to pure */
          /*  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           *    -->
           *  bc(r^{n-1}) dr^{n-1} */
          acc = nmove_restraineddown2pure(pi, pj, g, x, r2ij, gc->Zr[iens]);
          //printf("t %g, nmove- 2, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gc, iens - acc, pi, pj, t, "nmove_restraineddown2pure");
#endif
        } else if (gc->type[iens] == GCX_PURE) {
          /*  bc(R^n) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n */
          int ned = dg_randedge(g, &pi, &pj);

          if ( dg_connectedvs(g, mkbitsmask(g->n) ^ MKBIT(pi) ^ MKBIT(pj))
            && r2ij[pj][pi] < dblsqr(gc->rc[iens - 1] * gc->sr[iens - 1])
            && rnd0() < gc->Zr[iens] * ned ) {
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
                  1. / gc->sr[iens - 1], 1., 0);
          }
#ifdef CHECK
          //printf("acc %d, t %g, iens %d-->%d, type %d, pi %d, pj %d, rc %g, sr %g\n", acc, t, iens, iens - 1, gc->type[iens], pi, pj, gc->rc[iens - 1], gc->sr[iens - 1]);
          check(g, x, r2ij, gc, iens - acc, pi, pj, t, "nmove_purescaledown");
#endif
        } else { /* restrained down to a more restrained state */
          /*  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           * with Rc > rc */
          if (r2ij[pj][pi] < dblsqr(gc->rc[iens - 1] * gc->sr[iens - 1]))
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
                  1. / gc->sr[iens - 1], 1., 0);
          //printf("t %g, nmove- 1, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gc, iens - acc, pi, pj, t, "nmove_restraineddown2pure");
#endif
        }
        gc->ndown[iens][0] += 1;
        gc->ndown[iens][1] += acc;
        if ( acc ) iens--;
      }
    }
    else /* configuration sampling */
    {
      ctot += 1;
      if (iens % mtiers == 0) { /* pure state sampling */
        int i;
        //printf("t %g, cmove 0\n", t);
        BCSTEPR2(acc, i, g->n, g, ng, x, xi, r2ij, r2i, mcamp, gaussdisp);
        cacc += acc;
#ifdef CHECK
        check(g, x, r2ij, gc, iens, pi, pj, t, "cfg");
#endif
      } else { /* distance-restrained sampling */
        //printf("t %g, cmove 1\n", t);
        cacc += bcrstep(pi, pj, g, ng, x, xi, r2ij, r2i,
            gc->rc[iens], mcamp, gaussdisp);
        prtot += 1;
        pracc += changenapair(&pi, &pj, g, r2ij, gc->rc[iens]);
#ifdef CHECK
        check(g, x, r2ij, gc, iens, pi, pj, t, "cfgr");
#endif
      }
    }
STEP_END:
    /* center the structure */
    if (it % nstcom == 0) shiftr2ij(g, x, r2ij);
    gc->hist[iens] += 1;
    if (gc->type[iens] == GCX_PURE && rnd0() * nsted < 1)
      gc_accumdata(gc, g, t, nstcs, nstfb);
    if (ieql) { /* equilibration */
      if (it % nequil == 0) {
        gc_update(gc, mindata, updzr, gc->Zr, gc->rc, gc->sr);
        gc_computeZ(gc, gc->Zr, gc->rc, gc->sr);
        gc_saveZrr(gc, "Zrr.tmp", gc->Zr, gc->rc, gc->sr);
        gc_print(gc, 0, NULL, NULL, NULL);
        printf("equilibration stage %d/%d\n", ieql, neql);
        gc_cleardata(gc);
        t = 0; it = 0; /* reset time */
        if (++ieql >= neql) /* stop equilibration */
          ieql = 0;
      }
    } else { /* production */
      if (it % nstsave == 0 || t > nsteps - .5) {
        gc_update(gc, mindata, updzr, gc->Zr1, gc->rc1, gc->sr1);
        gc_computeZ(gc, gc->Zr1, gc->rc1, gc->sr1);
        if (restart) { /* don't write the new Zr for a restartable simulation */
          gc_saveZrr(gc, fnZrr, gc->Zr, gc->rc, gc->sr);
        } else {
          gc_saveZrr(gc, fnZrr, gc->Zr1, gc->rc1, gc->sr1);
        }
        gc_saveZr(gc, fnZr);
        it = 0;
      }
    }
  } /* main loop ends here */

  gc_print(gc, 0, gc->Zr1, gc->rc1, gc->sr1);
  printf("cacc %g, pracc %g\n", 1.*cacc/ctot, 1.*pracc/prtot);
  //{ int i; for (i = 0; i < g->n; i++) rvn_print(x[i], "x", "%+8.3f", 1); }
  gc_close(gc);
  dg_close(g);
  dg_close(ng);
  dg_close(g1);
  dg_close(ng1);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mcgcr(nmin, nmax, mtiers, nsteps, mcamp, neql, nequil, nstsave);
  mtsave(NULL);
  return 0;
}

