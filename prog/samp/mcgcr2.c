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

int updrc = 0; /* update rc */
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
typedef struct {
  int nmin, nmax;
  int m; /* number of ensembles for each n */
  int ens0; /* starting ensemble index */
  int nens; /* number of ensembles */
  int *type; /* 0 for pure */
  int *n; /* number of vertices */
  double *Zr; /* Z[iens] / Z[iens - 1] */
  double *Zr1; /* alternative buffer for Zr */
  real *rc, *rc2, *vol; /* volume of the insertion */
  real *sr; /* scaling of the distance between the first and last vertex */
  real *rc1, *sr1; /* alternative buffer for rc and sr */
  double *Z; /* partition function */
  double *hist; /* histogram */
  double (*ndown)[2]; /* number of upward (vertex adding) moves */
  double (*nup)[2]; /* number of downward (vertex removing) moves */
  double (*nedg)[2]; /* number of edges */
  double (*ncsp)[2]; /* no clique separator */
  double (*fbsm)[3]; /* hard-sphere weight */
  double (*ngpr)[2]; /* number of good pairs */
  double *B;
} gcx_t;



/* clear data */
static void gcx_cleardata(gcx_t *gcx)
{
  int i;

  for (i = 0; i < gcx->nens; i++) {
    gcx->hist[i] = 0;
    gcx->nup[i][0] = gcx->ndown[i][0] = 1e-20;
    gcx->nup[i][1] = gcx->ndown[i][1] = 0;
  }
  for (i = 0; i <= gcx->nmax; i++) {
    gcx->nedg[i][0] = 1e-20;
    gcx->nedg[i][1] = 0;
    gcx->ncsp[i][0] = 1e-20;
    gcx->ncsp[i][1] = 0;
    gcx->fbsm[i][0] = 1e-20;
    gcx->fbsm[i][1] = gcx->fbsm[i][2] = 0;
    gcx->ngpr[i][0] = 1e-20;
    gcx->ngpr[i][1] = 0;
  }
}



#define GCX_PURE        0
#define GCX_RESTRAINED  1

INLINE gcx_t *gcx_open(int nmin, int nmax, int m,
    real rc0, real sr0)
{
  gcx_t *gcx;
  int i, n;

  xnew(gcx, 1);
  gcx->nmin = nmin;
  gcx->nmax = nmax;
  gcx->m = m;
  gcx->nens = 1 + nmax * m;
  gcx->ens0 = nmin * m;
  xnew(gcx->type, gcx->nens);
  xnew(gcx->n, gcx->nens);
  xnew(gcx->rc, gcx->nens);
  xnew(gcx->rc2, gcx->nens);
  xnew(gcx->vol, gcx->nens);
  xnew(gcx->sr, gcx->nens);
  /* temporary setting */
  if (m > 2 && rc0 >= 1) rc0 = 0.8;
  for (i = 0; i < gcx->nens; i++) {
    gcx->n[i] = (i + m - 1) / m;
    /* type == 0 mean GCX_PURE */
    gcx->type[i] = (i % m);
    //printf("init. i %d, n %d, rc %g\n", i, gcx->n[i], gcx->rc[i]);
  }
  for (i = 0; i < gcx->nens; i++) {
    n = gcx->n[i];
    if ( gcx->type[i] == GCX_PURE ) {
      gcx->rc[i] = 1;
    } else {
      /* gcx->rc[i] = pow(rc0, 1.*(m - i % m)/(m - 1)); */
      if (m <= 2) { /* a single intermediate state */
        gcx->rc[i] = rc0;
      } else if (m == 3 || n <= 6) { /* two intermediate states */
        /* linear interpolation, excluding the value 1 */
        gcx->rc[i] = rc0 + (1. - rc0) * (i % m - 1) / (m - 1);
      } else if (m >= 4) {
        /* linear interpolation, including the boundary value 1 */
        gcx->rc[i] = rc0 + (1. - rc0) * (i % m - 1) / (m - 2);
      }
    }
    gcx->rc2[i] = gcx->rc[i] * gcx->rc[i];
    gcx->vol[i] = pow(gcx->rc[i], D);
  }

  for (i = 0; i < gcx->nens; i++) {
    if ( i == gcx->nens - 1 || gcx->type[i] == GCX_PURE ) {
      gcx->sr[i] = 1;
    } else if ( sr0 > 1 && i > 0 && gcx->type[i - 1] == GCX_PURE ) {
      /* the given value is only used for the first rc */
      gcx->sr[i] = sr0;
    } else { /* guess from the ratios of the successive radii */
      real den = gcx->rc[i];
      if (den > 1) gcx->sr[i] = 1;
      else gcx->sr[i] = gcx->rc[i + 1] / den;
    }
  }
  xnew(gcx->Zr, gcx->nens);
  xnew(gcx->Z, gcx->nens);
  for (i = 0; i < gcx->nens; i++) {
    if (gcx->type[i] == GCX_PURE)
      gcx->Zr[i] = 1./gcx->n[i]; /* inverse number of nonarticulated pairs */
    else gcx->Zr[i] = 1;
    gcx->Z[i] = 1;
  }

  /* alternative buffer */
  xnew(gcx->Zr1, gcx->nens);
  xnew(gcx->rc1, gcx->nens);
  xnew(gcx->sr1, gcx->nens);

  xnew(gcx->hist, gcx->nens);
  xnew(gcx->ndown, gcx->nens);
  xnew(gcx->nup, gcx->nens);

  /* essential statistics */
  xnew(gcx->nedg, gcx->nmax + 1);
  xnew(gcx->ncsp, gcx->nmax + 1);
  xnew(gcx->fbsm, gcx->nmax + 1);
  xnew(gcx->ngpr, gcx->nmax + 1);
  xnew(gcx->B, gcx->nmax + 1);
  for (i = 0; i <= gcx->nmax; i++)
    gcx->B[i] = 1;
  gcx_cleardata(gcx);
  return gcx;
}



INLINE void gcx_close(gcx_t *gcx)
{
  free(gcx->type);
  free(gcx->n);
  free(gcx->rc);
  free(gcx->rc2);
  free(gcx->vol);
  free(gcx->sr);
  free(gcx->Zr);
  free(gcx->Z);
  free(gcx->Zr1);
  free(gcx->rc1);
  free(gcx->sr1);
  free(gcx->hist);
  free(gcx->ndown);
  free(gcx->nup);
  free(gcx->nedg);
  free(gcx->ncsp);
  free(gcx->fbsm);
  free(gcx->ngpr);
  free(gcx->B);
  free(gcx);
}



/* accumulate data */
static void gcx_accumdata(gcx_t *gcx, const dg_t *g, double t,
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
  gcx->nedg[n][0] += 1;
  gcx->nedg[n][1] += ned;

  if (n <= nlookup || rnd0() * nstcs / nsted < 1) {
    /* check if the graph has a clique separator */
    if (n <= nlookup) { /* lookup table */
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
    gcx->ncsp[n][0] += 1;
    gcx->ncsp[n][1] += ncs;

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
      gcx->fbsm[n][0] += 1;
      gcx->fbsm[n][1] += fb;
      gcx->fbsm[n][2] += abs(fb);
    }
  }
}



/* update parameters */
static void gcx_update(gcx_t *gcx, double mindata, int updrc,
    double *Zr, real *rc, real *sr)
{
  int i;
  double r;

  if (Zr != gcx->Zr) memcpy(Zr, gcx->Zr, sizeof(gcx->Zr[0]) * gcx->nens);
  if (rc != gcx->rc) memcpy(rc, gcx->rc, sizeof(gcx->rc[0]) * gcx->nens);
  if (sr != gcx->sr) memcpy(sr, gcx->sr, sizeof(gcx->sr[0]) * gcx->nens);
  for (i = 1; i < gcx->nens; i++) {
    /* make sure enough data points and not to change the exact data */
    if (gcx->ndown[i][1] >= mindata && gcx->nup[i-1][1] >= mindata) {
      r = (gcx->nup[i-1][1] / gcx->nup[i-1][0])
        / (gcx->ndown[i][1] / gcx->ndown[i][0]);
      if (gcx->type[i - 1] == GCX_PURE) {
        /* a pure state and an upper distance-restrained state */
        if (updrc) { /* updating rc is disabled by default */
          /* NOTE: we update rc such that the forward acceptance ratio equals
           * the backward one, but it may settle in a sub-optimal value
           * for the forward (adding vertex) acceptance ratio is not monotonic
           * with respect to the magnitude of rc */
          /* it appear r^(1/4) converges better than r^(1/D), but 1./D
           * is more stable, and we don't want to change rc much */
          rc[i] *= pow(r, 1./D);
        } else {
          Zr[i] *= r;
        }
      } else if (gcx->type[i] == GCX_PURE) {
        /* we choose Zr as 1/|ngpr| to increase acceptance */
        int n = gcx->n[i];
        double nZr = Zr[i];
        //if (gcx->ngpr[n][1] > mindata) nZr = gcx->ngpr[n][0] / gcx->ngpr[n][1];
        sr[i - 1] *= pow(r * Zr[i]/nZr, 1./D);
        //printf("n %d, updating Zr[%d] from %g to %g\n", n, i, Zr[i], nZr);
        Zr[i] = nZr;
      } else {
        sr[i - 1] *= pow(r, 1./D);
      }
      //printf("update i %3d, n %2d, r %.7f, rc %.7f, sr %.7f\n", i, gcx->n[i], r, rc[i], sr[i - 1]);
    }
  }
}



/* compute the partition function and virial coefficients */
static void gcx_computeZ(gcx_t *gcx,
    const double *Zr, const real *rc, const real *sr)
{
  double fac = 1, x = 1, fbav;
  int n, i;

  for (i = 1; i < gcx->nens; i++) {
    n = gcx->n[i];
    /* compute the partition function */
    x *= Zr[i];
    if (gcx->type[i - 1] != GCX_PURE) {
      x *= pow(sr[i - 1], D);
    } else {
      x *= pow(rc[i], D);
    }
    //printf("compute Z, i %d, n %d, rc %.7f, sr %.7f, Z %.6e\n", i, n, rc[i], sr[i - 1], x);
    /* use the the exact result if available */
    if (gcx->type[i] == GCX_PURE) {
      x *= n * (n - 1) * .5; /* pair generating bias */
      if (n == 4 && x <= 12) {
        printf("Z4 %g vs %g(ref), Zri %g, (sr*rc)^D %g\n", x, Z4rat(D), gcx->Zr[i], pow(gcx->sr[i-1] * gcx->rc[i-1], D));
      }
      if (n <= 2) x = 1;
      else if (n == 3) x = Z3rat(D);
      else if (n == 4 && D <= 12) x = Z4rat(D);
    }
    gcx->Z[i] = x;
    if (gcx->type[i] == GCX_PURE) {
      if (n >= 2) fac *= 2./n;
      fbav = gcx->fbsm[n][1] / gcx->fbsm[n][0];
      /* Bn/B2^(n-1) = (-n) 2^(n - 1) / n! Zn/Z2^(n-1) < fb > */
      gcx->B[n] = (1. - n) * fac * x * fbav;
      //printf("B[%d] = %g\n", n, gcx->B[n]);
    }
  }
}



/* print out a summary */
static void gcx_print(const gcx_t *gcx, int compact,
    double *Zr, double *rc, double *sr)
{
  int i;

  if (Zr == NULL) Zr = gcx->Zr;
  if (rc == NULL) rc = gcx->rc;
  if (sr == NULL) sr = gcx->sr;
  for (i = gcx->ens0; i < gcx->nens; i++) {
    printf("%3d %2d %9.6f %8.6f %8.6f ",
       i, gcx->n[i], Zr[i], rc[i], sr[i]);
    if ( !compact ) {
      printf("%12.0f %.5f %.5f ", gcx->hist[i],
        gcx->nup[i-1][1] / gcx->nup[i-1][0],
        gcx->ndown[i][1] / gcx->ndown[i][0]);
      if (gcx->type[i] == GCX_PURE) {
        int n = gcx->n[i];
        printf("%6.3f %6.3f %8.6f %+9.7f %+.5e",
          gcx->ngpr[n][1] / gcx->ngpr[n][0],
          gcx->nedg[n][1] / gcx->nedg[n][0],
          gcx->ncsp[n][1] / gcx->ncsp[n][0],
          gcx->fbsm[n][1] / gcx->fbsm[n][0],
          gcx->B[n]);
      }
    }
    printf("\n");
  }
}



#define mkfnZrdef(fn, fndef, d, m, n) \
  if (fn == NULL) { sprintf(fndef, "ZrD%dr%dn%d.dat", d, m, n); fn = fndef; }



/* save a compact (pure-state) Zr file, mimic the output of mcgc.c */
static int gcx_saveZr(gcx_t *gcx, const char *fn)
{
  FILE *fp;
  int i, n;
  char fndef[64];
  double Zr = 1, prevZ = 1, tot = 0;

  mkfnZrdef(fn, fndef, D, gcx->m - 1, gcx->nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d V2 %d\n", D, gcx->nmax, nedxmax);
  for (i = 1; i < gcx->nens; i++) {
    if (gcx->type[i] != GCX_PURE) continue;
    n = gcx->n[i];
    /* start with the third virial coefficient */
    if (n <= 2) continue;
    Zr = gcx->Z[i] / prevZ;
    fprintf(fp, "%3d %18.14f %20.14e %14.0f %14.0f %14.0f %.14f %.14f "
        "%14.0f %14.0f %14.0f %18.14f %16.14f %16.14f %+17.14f %+20.14e\n",
        n, Zr, gcx->Z[i], gcx->hist[i],
        gcx->nup[i-1][0], gcx->ndown[i][0],
        gcx->nup[i-1][1] / gcx->nup[i-1][0],
        gcx->ndown[i][1] / gcx->ndown[i][0],
        gcx->nedg[n][0], gcx->ncsp[n][0], gcx->fbsm[n][0],
        gcx->nedg[n][1] / gcx->nedg[n][0],
        gcx->ncsp[n][1] / gcx->ncsp[n][0],
        gcx->fbsm[n][2] / gcx->fbsm[n][0],
        gcx->fbsm[n][1] / gcx->fbsm[n][0],
        gcx->B[n]);
    tot += gcx->hist[i];
    prevZ = gcx->Z[i];
  }
  fclose(fp);
  printf("saved compact Zr to %s, total %g\n", fn, tot);
  return 0;
}



#define mkfnZrrdef(fn, fndef, d, m, n) \
  if (fn == NULL) { sprintf(fndef, "ZrrD%dr%dn%d.dat", d, m, n); fn = fndef; }



/* save all data to file */
static int gcx_saveZrr(gcx_t *gcx, const char *fn,
    const double *Zr, const double *rc, const double *sr)
{
  FILE *fp;
  int i, n;
  char fndef[64];
  double tot = 0;

  mkfnZrrdef(fn, fndef, D, gcx->m - 1, gcx->nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#R %d %d %d %d %d %d V2 %d\n",
      D, gcx->nens, gcx->ens0, gcx->nmin, gcx->nmax, gcx->m, nedxmax);
  for (i = 1; i < gcx->nens; i++) {
    fprintf(fp, "%4d %d %3d %20.14e %20.14e "
        "%18.14f %18.14f %14.0f %14.0f %14.0f %.14f %.14f ",
        i, gcx->type[i], gcx->n[i], Zr[i], gcx->Z[i], rc[i], sr[i],
        gcx->hist[i],
        gcx->nup[i-1][0], gcx->ndown[i][0],
        gcx->nup[i-1][1] / gcx->nup[i-1][0],
        gcx->ndown[i][1] / gcx->ndown[i][0]);
    if (gcx->type[i] == GCX_PURE) {
      n = gcx->n[i];
      fprintf(fp, "%14.0f %14.0f %14.0f ",
        gcx->nedg[n][0], gcx->ncsp[n][0], gcx->fbsm[n][0]);
      fprintf(fp, "%18.14f %16.14f %16.14f %+17.14f %+20.14e",
        gcx->nedg[n][1] / gcx->nedg[n][0],
        gcx->ncsp[n][1] / gcx->ncsp[n][0],
        gcx->fbsm[n][2] / gcx->fbsm[n][0],
        gcx->fbsm[n][1] / gcx->fbsm[n][0],
        gcx->B[n]);
    }
    fprintf(fp, "\n");
    tot += gcx->hist[i];
  }
  fclose(fp);
  printf("saved Zrr to %s, tot %g\n", fn, tot);
  return 0;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int gcx_load(gcx_t *gcx, const char *fn, int loadall)
{
  FILE *fp;
  int i = -1, nens, d = 0, n, tp, ens0 = 0, n0, m = 0, ver = 0, next;
  char s[512], fndef[64], sver[16], *p;

  mkfnZrrdef(fn, fndef, D, gcx->m - 1, gcx->nmax);
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
    || d != D || m != gcx->m || ens0 != gcx->ens0 ) {
    fprintf(stderr, "%s: D %d vs. %d, m %d vs %d, ens0 %d vs %d\n%s",
        fn, d, D, m, gcx->m, ens0, gcx->ens0, s);
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
        || i != d || tp != gcx->type[i] || n != gcx->n[i]) {
      fprintf(stderr, "%s ends on line %d\n%s", fn, i, s);
      break;
    }
    if (Zr >= 0) gcx->Zr[i] = Zr;
    if (rc >= 0) {
      gcx->rc[i] = rc;
      gcx->rc2[i] = rc * rc;
      gcx->vol[i] = pow(rc, D);
    }
    if (sr >= 0) {
      if (ver >= 2) gcx->sr[i] = sr;
      else gcx->sr[i - 1] = sr;
    }

    if (loadall) { /* try to get additional data */
      p = s + next;
      if (5 != sscanf(p, "%lf%lf%lf%lf%lf%n",
            &gcx->hist[i], &gcx->nup[i-1][0], &gcx->ndown[i][0],
            &gcx->nup[i-1][1], &gcx->ndown[i][1], &next) ) {
        fprintf(stderr, "%s corrupted on line %d\n%s", fn, i, s);
        break;
      }
      if (gcx->nup[i-1][0] <= 0) gcx->nup[i-1][0] = 1e-20;
      if (gcx->ndown[i][0] <= 0) gcx->ndown[i][0] = 1e-20;
      gcx->nup[i-1][1] *= gcx->nup[i-1][0];
      gcx->ndown[i][1] *= gcx->ndown[i][0];
      if (gcx->type[i] == GCX_PURE) {
        p += next;
        if (7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf",
              &gcx->nedg[n][0], &gcx->ncsp[n][0], &gcx->fbsm[n][0],
              &gcx->nedg[n][1], &gcx->ncsp[n][1],
              &gcx->fbsm[n][2], &gcx->fbsm[n][1]) ) {
          fprintf(stderr, "%s corrupted on line %d\n%s", fn, i, s);
          break;
        }
        if (gcx->nedg[n][0] <= 0) gcx->nedg[n][0] = 1e-20;
        if (gcx->ncsp[n][0] <= 0) gcx->ncsp[n][0] = 1e-20;
        if (gcx->fbsm[n][0] <= 0) gcx->fbsm[n][0] = 1e-20;
        gcx->nedg[n][1] *= gcx->nedg[n][0];
        gcx->ncsp[n][1] *= gcx->ncsp[n][0];
        gcx->fbsm[n][1] *= gcx->fbsm[n][0];
        gcx->fbsm[n][2] *= gcx->fbsm[n][0];
      }
    }
  }
  fclose(fp);
  if (i >= nens) /* recompute the partition function */
    gcx_computeZ(gcx, gcx->Zr, gcx->rc, gcx->sr);
  printf("loaded Zr from %s, version %d\n", fn, ver);
  return i;
}



/* compute the distance matrix */
INLINE void calcr2ij(real r2ij[][DG_NMAX], rvn_t *x, int n)
{
  int i, j;

  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++)
      r2ij[j][i] = r2ij[i][j] = rvn_dist2(x[i], x[j]);
}



/* build graph from the distance matrix r2ij (lower block) */
INLINE void mkgraphr2ij(dg_t *g, real r2ij[][DG_NMAX], real s, int n)
{
  int i, j;
  real s2 = s * s;

  g->n = n;
  dg_empty(g);
  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++)
      if (r2ij[i][j] < s2) {
        DG_LINK(g, i, j);
      }
}



/* construct a new graph with a single vertex i displaced
 * and save the corresponding array of displacements */
#define UPDGRAPHR2(i, nv, g, ng, x, xi, r2, r2i) { int j; \
  ng->n = (nv); \
  dg_copy(ng, g); \
  for ( (r2i)[i] = 0, j = 0; j < (nv); j++) { \
    if ( j == i ) continue; \
    if ( ((r2i)[j] = rvn_dist2(xi, x[j])) < (r2) ) { \
      DG_LINK(ng, i, j); \
    } else { \
      DG_UNLINK(ng, i, j); \
    } \
  } }



/* update the matrix r2ij */
#define UPDR2IJ(r2ij, r2i, nv, i, j) { \
  for (j = 0; j < (nv); j++) { \
    if (j != i) (r2ij)[i][j] = (r2ij)[j][i] = (r2i)[j]; \
  } }



/* single-vertex move in a pure state */
INLINE int bcstep(dg_t *g, dg_t *ng, rvn_t *x, real *xi,
    real r2ij[][DG_NMAX], real *r2i, real amp, int gauss)
{
  int i, j, n = g->n;

  DISPRNDI(i, n, x, xi, amp, gauss);
  UPDGRAPHR2(i, n, g, ng, x, xi, 1, r2i);
  if ( dg_biconnected(ng) ) {
    rvn_copy(x[i], xi);
    dg_copy(g, ng);
    UPDR2IJ(r2ij, r2i, n, i, j);
    return 1;
  }
  return 0;
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
    UPDR2IJ(r2ij, r2i, n, i, j);
    return 1;
  }
  return 0;
}



/* attach a vertex to a random vertex i0, become a restrained state
 *  pure: bc(r^n)
 *  restrained: step(Rc - R_{i0, n}) bc(R^{n+1}) c(R^{n+1}\{i0, n})
 *  transition state:
 *    bc(r^n) step(Rc - R_{i0, n}) bi(R^{n+1}) d R^{n+1}
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
 *  restrained: bc(r^n) step(rc - r_{i0, j0}) c(r^{n}\{i0,j0})
 *  pure: bc(r^{n-1})
 *  transition state:  bc(r^{n-1}) step(rc - r_{i0, j0}) bi(r^n)
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



static int napr_[DG_NMAX*(DG_NMAX - 1)/2];

/* list all non-articulated pairs within rc */
INLINE int listpairs(dg_t *g, real r2ij[][DG_NMAX], real rc2)
{
  int npr = 0, i, j, n = g->n;
  code_t mask = mkbitsmask(n);

  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++)
      if ( r2ij[i][j] < rc2
        && dg_connectedvs(g, mask ^ MKBIT(i) ^ MKBIT(j)) )
        napr_[npr++] = i * DG_NMAX + j;
  return npr;
}



/* randomly choose a non-articulated pair within rc
 * return the number of qualified pairs */
INLINE int nmove_choosepair(int *pi, int *pj, dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int npr = 0, ipr, id;

  npr = listpairs(g, r2ij, rc * rc);
  if (npr && pi != NULL && pj != NULL) {
    id = napr_[ (int) (rnd0() * npr) ];
    *pi = id / DG_NMAX;
    *pj = id - (*pi) * DG_NMAX;
  }
  return npr;
}



#define changepair(i, j, g, r2ij, r2) changepair_heatbath(i, j, g, r2ij, r2)

/* change the distance-restrained pair */
INLINE int changepair_metro(int *i, int *j, dg_t *g,
    real r2ij[][DG_NMAX], real rc2)
{
  int ni, nj, found = 0;

  ni = randpair(g->n, &nj);
  if (ni == *i && nj == *j) return 0;
  if (ni == *j && nj == *i) {
    *i = ni;
    *j = nj;
    return 0;
  }
  if ( dg_connectedvs(g, mkbitsmask(g->n) ^ MKBIT(ni) ^ MKBIT(nj))
    && r2ij[ni][nj] < rc2 ) {
    *i = ni;
    *j = nj;
    return 1;
  }
  return found;
}



/* change the distance-restrained pair, heat bath algorithm */
INLINE int changepair_heatbath(int *pi, int *pj, dg_t *g,
    real r2ij[][DG_NMAX], real rc2)
{
  int npr, ipr, ipr0, n = g->n, pid;

  npr = listpairs(g, r2ij, rc2);
  if (npr <= 1) { /* we can at least swap i and j */
    pid = *pi, *pi = *pj, *pj = pid;
    return 0;
  }
  pid = (*pi > *pj) ? (*pi * DG_NMAX + *pj) : (*pj * DG_NMAX + *pi);
  for (ipr0 = 0; ipr0 < npr; ipr0++)
    if (napr_[ipr0] == pid) break;
  die_if (ipr0 >= npr, "cannot find pair %d, %d, npr %d\n", *pi, *pj, npr);
  ipr = (int) (rnd0() * (npr - 1));
  if (ipr >= ipr0) ipr++;
  pid = napr_[ipr];
  *pi = pid / DG_NMAX;
  *pj = pid - (*pi) * DG_NMAX;
  if (rnd0() < 0.5) {
    pid = *pi, *pi = *pj, *pj = pid;
  }
  return 1;
}



/* scale the distance between x[i0] and x[j0] and test the biconnectivity
 * used for upwards moves to a pure state */
INLINE int nmove_scaleup2pure(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xj0, real r2ij[][DG_NMAX], real *r2i,
    real rc, real sr, real Zr)
{
  int k, n = g->n;

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

  /* count the number of non-articulated pairs */
  k = nmove_choosepair(NULL, NULL, ng, r2ij, sr * rc);
  Zr *= 1. / k;
  if (Zr < 1 && rnd0() >= Zr) return 0;

  dg_copy(g, ng);
  rvn_copy(x[j0], xj0);
  UPDR2IJ(r2ij, r2i, n, j0, k);
  return 1;
}



/* scale the distance between x[i0] and x[j0] and test the biconnectivity */
INLINE int nmove_scale(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xj0, real r2ij[][DG_NMAX], real *r2i,
    real sr, real Zr)
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
  if ( dg_biconnected(ng) ) {
    dg_copy(g, ng);
    rvn_copy(x[j0], xj0);
    UPDR2IJ(r2ij, r2i, n, j0, k);
    return 1;
  }
  return 0;
}



/* center the configuration, recompute r2ij for precision */
INLINE void center(const dg_t *g, rvn_t *x, real r2ij[][DG_NMAX])
{
  int i, j, n = g->n;

  /* put the first vertex to the center */
  for (i = 1; i < n; i++) rvn_dec(x[i], x[0]);
  rvn_zero(x[0]);
  /* recompute the matrix r2ij */
  for (i = 1; i < n; i++)
    for (j = 0; j < i; j++) {
      real r2 = rvn_dist2(x[i], x[j]);
      int lnk = dg_linked(g, i, j);
      die_if ( fabs(r2 - r2ij[i][j]) > 1e-6
            || (r2 < 1 && !lnk) || (r2 >= 1 && lnk),
        "r2ij[%d][%d] = %g vs. %g, link %d\n",
        i, j, r2, r2ij[i][j], lnk);
      r2ij[i][j] = r2ij[j][i] = r2;
    }
}



#ifdef CHECK
/* check if everything is okay */
INLINE int check(const dg_t *g, rvn_t *x, real (*r2ij)[DG_NMAX],
    const gcx_t *gcx, int iens, int i0, int j0,
    double t, const char *msg)
{
  int i, j, n = g->n, err = 0, type = gcx->type[iens];
  real r2, rc2 = gcx->rc2[iens];

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
          i0, j0, n, r2ij[i0][j0], rc2, gcx->rc[iens]);
      err++;
    }
    if ( !dg_connectedvs(g, mkbitsmask(g->n) ^ MKBIT(i0) ^ MKBIT(j0)) ) {
      fprintf(stderr, "the pair (%d, %d) is articulated\n", i0, j0);
      err++;
    }
  }
  die_if (err, "t %g, iens %d, type %d, n %d, message %s, err %d\n",
      t, iens, gcx->type[iens], gcx->n[i], msg, err);
  return err;
}
#endif



static void mcgcr(int nmin, int nmax, int mtiers, double nsteps,
    real mcamp, int neql, int nequil, int nstsave)
{
  double t, ctot = 0, cacc = 0, pracc = 0;
  int iens, ensmax, ensmin, acc, it, ieql;
  int pi = 0, pj = 1; /* indices of the distance-restrained pair */
  gcx_t *gcx;
  dg_t *g, *ng, *g1, *ng1;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  real r2ij[DG_NMAX][DG_NMAX] = {{0}}, r2i[DG_NMAX] = {0};

  g = dg_open(nmax);
  ng = dg_open(nmax);
  g1 = dg_open(nmax);
  ng1 = dg_open(nmax);
  gcx = gcx_open(nmin, nmax, mtiers, rc0, sr0);
  gcx_load(gcx, fninp, restart && !bsim0);
  gcx_print(gcx, 1, NULL, NULL, NULL);
  /* initial n is nmin */
  calcr2ij(r2ij, x, nmin);

  ensmin = nmin * mtiers;
  ensmax = nmax * mtiers;
  iens = ensmin;
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
        if (gcx->type[iens] == GCX_PURE) { /* pure up to restrained */
          /*  bc(r^n) dr^n
           *    -->
           *  bc(R^{n+1}) step(Rc - R_{i0, n}) c(r^{n+1}\{i0, n}) dR^{n+1} */
          pj = g->n;
          acc = nmove_pureup2restrained(&pi, g, x, r2ij,
              gcx->rc[iens + 1], 1. / gcx->Zr[iens + 1]);
          //printf("t %g, n+ move 1, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gcx, iens + acc, pi, pj, t, "nmove_pureup2restrained");
#endif
          die_if (g->n < nmin || g->n > nmax, "p2m BAD n %d, t %g, iens %d\n", g->n, t, iens);
        } else if (gcx->type[iens + 1] == GCX_PURE) { /* restrained up to pure */
          /*  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           *    -->
           *  bc(R^n) dR^n */
          acc = nmove_scaleup2pure(pi, pj, g, ng, x, xi, r2ij, r2i,
            gcx->rc[iens], gcx->sr[iens], 1. / gcx->Zr[iens + 1]);
        } else { /* restrained up to a restrained state */
          /*  c(r^n\{i0, j0}) step(rc - r_{i0, j0}) bc(r^n) dr^n
           *    -->
           *  c(R^n\{i0, j0}) step(Rc - R_{i0, j0}) bc(R^n) dR^n
           * where rc < Rc */
          if ( r2ij[pi][pj] * dblsqr(gcx->sr[iens]) < gcx->rc2[iens + 1] )
            /* move after a scaling between the first and last vertices */
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
              gcx->sr[iens], 1. / gcx->Zr[iens + 1]);
          //printf("t %g, n+ move 2, acc %d, iens %d, n %d, n1 %d\n", t, acc, iens, g->n, g1->n);
#ifdef CHECK
          check(g, x, r2ij, gcx, iens + acc, pi, pj, t, "nmove_scaleup");
#endif
        }
        gcx->nup[iens][0] += 1;
        gcx->nup[iens][1] += acc;
        if ( acc ) iens++;
        die_if (g->n < nmin || g->n > nmax, "BAD n %d, t %g, iens %d\n", g->n, t, iens);
      }
      else /* decrease the ensemble index */
      {
        if (iens <= ensmin) goto STEP_END;

        acc = 0;
        if (gcx->type[iens - 1] == GCX_PURE) { /* restrained down to pure */
          /*  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           *    -->
           *  bc(r^{n-1}) dr^{n-1} */
          acc = nmove_restraineddown2pure(pi, pj, g, x, r2ij,
              gcx->Zr[iens]);
          //printf("t %g, nmove- 2, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gcx, iens - acc, pi, pj, t, "nmove_restraineddown2pure");
#endif
        } else if (gcx->type[iens] == GCX_PURE) {
          /*  bc(R^n) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n */
          int ngpr = nmove_choosepair(&pi, &pj, g, r2ij,
              gcx->rc[iens - 1] * gcx->sr[iens - 1]);
          int n = gcx->n[iens];
          if (ngpr) { /* ngpr is the number of non-articulated pairs */
            /* move after a scaling of two vertices */
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
                1. / gcx->sr[iens - 1], gcx->Zr[iens] * ngpr);
            gcx->ngpr[n][0] += 1;
            gcx->ngpr[n][1] += ngpr;
          }
          //printf("t %g, nmove- 1, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          //printf("acc %d, t %g, iens %d-->%d, type %d, pi %d, pj %d, rc %g, sr %g\n", acc, t, iens, iens - 1, gcx->type[iens], pi, pj, gcx->rc[iens - 1], gcx->sr[iens - 1]);
          check(g, x, r2ij, gcx, iens - acc, pi, pj, t, "nmove_purescaledown");
#endif
        } else { /* restrained down to restrained */
          /*  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           * with Rc > rc */
          if (r2ij[pj][pi] < dblsqr(gcx->rc[iens - 1] * gcx->sr[iens - 1]))
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i,
                  1. / gcx->sr[iens - 1], gcx->Zr[iens]);
          //printf("t %g, nmove- 1, acc %d, iens %d, n %d\n", t, acc, iens, g->n);
#ifdef CHECK
          check(g, x, r2ij, gcx, iens - acc, pi, pj, t, "nmove_restrictedscaledown");
#endif
        }
        gcx->ndown[iens][0] += 1;
        gcx->ndown[iens][1] += acc;
        if ( acc ) iens--;
      }
    }
    else /* configuration sampling */
    {
      ctot += 1;
      if (iens % mtiers == 0) { /* pure state sampling */
        //printf("t %g, cmove 0\n", t);
        cacc += bcstep(g, ng, x, xi, r2ij, r2i, mcamp, gaussdisp);
#ifdef CHECK
        check(g, x, r2ij, gcx, iens, pi, pj, t, "cfg");
#endif
      } else { /* distance-restrained sampling */
        //printf("t %g, cmove 1\n", t);
        cacc += bcrstep(pi, pj, g, ng, x, xi, r2ij, r2i,
            gcx->rc[iens], mcamp, gaussdisp);
        pracc += changepair(&pi, &pj, g, r2ij, gcx->rc2[iens]);
#ifdef CHECK
        check(g, x, r2ij, gcx, iens, pi, pj, t, "cfgr");
#endif
      }
    }
STEP_END:
    /* center the structure */
    gcx->hist[iens] += 1;
    if (gcx->type[iens] == GCX_PURE && rnd0() * nsted < 1)
      gcx_accumdata(gcx, g, t, nstcs, nstfb);
    if (ieql) { /* equilibration */
      if (it % nequil == 0) {
        gcx_update(gcx, mindata, updrc, gcx->Zr, gcx->rc, gcx->sr);
        gcx_computeZ(gcx, gcx->Zr, gcx->rc, gcx->sr);
        gcx_print(gcx, 0, NULL, NULL, NULL);
        printf("equilibration stage %d/%d\n", ieql, neql);
        gcx_cleardata(gcx);
        t = 0; it = 0; /* reset time */
        if (++ieql >= neql) /* stop equilibration */
          ieql = 0;
      }
    } else { /* production */
      if (it % nstcom == 0)
        center(g, x, r2ij);
      if (it % nstsave == 0 || t > nsteps - .5) {
        gcx_update(gcx, mindata, updrc, gcx->Zr1, gcx->rc1, gcx->sr1);
        gcx_computeZ(gcx, gcx->Zr1, gcx->rc1, gcx->sr1);
        if (restart) { /* don't write the new Zr for a restartable simulation */
          gcx_saveZrr(gcx, fnZrr, gcx->Zr, gcx->rc, gcx->sr);
        } else {
          gcx_saveZrr(gcx, fnZrr, gcx->Zr1, gcx->rc1, gcx->sr1);
        }
        gcx_saveZr(gcx, fnZr);
        it = 0;
      }
    }
  } /* main loop ends here */

  gcx_update(gcx, mindata, updrc, gcx->Zr, gcx->rc, gcx->sr);
  gcx_computeZ(gcx, gcx->Zr, gcx->rc, gcx->sr);
  gcx_print(gcx, 0, NULL, NULL, NULL);
  printf("cacc %g, pracc %g\n", 1.*cacc/ctot, 1.*pracc/ctot);
  //{ int i; for (i = 0; i < g->n; i++) rvn_print(x[i], "x", "%+8.3f", 1); }
  gcx_close(gcx);
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

