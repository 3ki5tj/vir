/* grand-canonical simulation to compute the partition function
 * of all biconnected diagrams
 * simple Widom insertion with Metropolis or heat-bath removal */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dgrjw.h"
#include "dgring.h"
#include "mcutil.h"

#ifdef MPI
#include "mpi.h"
#define MPI_MYREAL ( (sizeof(real) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT )
/* global variables */
MPI_Comm comm = MPI_COMM_WORLD;
#endif
#define MASTER 0
int inode = MASTER, nnodes = 1;



int nmin = 1; /* the minimal order of virial coefficients */
int nmax = DG_NMAX; /* the maximal order of virial coefficients */
real mcamp = 1.5f;
double nequil = 10000000;
double nsteps = 10000000;
double ratn = 0.5; /* frequency of n-moves */
  /* this version used a simple n-move, which is on average faster than
   * configurational sampling, so we should use a larger ratn */
real rc0 = 0; /* initial rc */
double mindata = 100; /* minimal # of data points to make update */

char *fninp = NULL; /* input file */
char *fnout = NULL; /* output file */
char *fnZrtmp = "Zrh.tmp"; /* temporary output file */
char *fnBring = NULL; /* ring contribution file */

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
int gaussdisp = 0;

int neql = -1; /* number of equilibration stages
                 == 0: disable
                 >  0: Zr, sr are updated after each round
                 <  0: Zr, sr are not updated, one round */
int nstsave = 1000000000; /* interval of saving data */

int restart = 0; /* restartable simulation */
int bsim0 = 0; /* first simulation */
int updrc = 1; /* update rc */

/* nmove type, 0: Metropolis, 1: heat-bath */
int nmvtype = 1;



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int norc = 0;
  double nstsav = nstsave;

  ao = argopt_open(0);
  argopt_add(ao, "-z",    "%d",   &nmin,    "minimal order of virial coefficients");
  argopt_add(ao, "-n",    "%d",   &nmax,    "maximal order of virial coefficients");
  argopt_add(ao, "-1",    "%lf",  &nsteps,  "number of simulation steps");
  argopt_add(ao, "-a",    "%r",   &mcamp,   "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-r",    "%lf",  &ratn,    "rate of particle moves");
  argopt_add(ao, "-c",    "%r",   &rc0,     "initial rc");
  argopt_add(ao, "-i",    NULL,   &fninp,   "input file");
  argopt_add(ao, "-o",    NULL,   &fnout,   "output file");
  argopt_add(ao, "-T",    NULL,   &fnZrtmp, "temporary output file");
  argopt_add(ao, "-I",    NULL,   &fnBring, "ring contribution file");
  argopt_add(ao, "-G",    "%d",   &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-A",    "%d",   &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-F",    "%d",   &nstfb,   "interval of computing fb");
  argopt_add(ao, "-X",    "%d",   &nedxmax, "maximal number of extra edges for computing fb");
  argopt_add(ao, "-Q",    "%lf",  &mindata, "minimal number of data points to warrant parameter update");
  argopt_add(ao, "-q",    "%lf",  &nstsav,  "interval of saving data");
  argopt_add(ao, "-E",    "%d",   &neql,    "number of equilibration rounds");
  argopt_add(ao, "-0",    "%lf",  &nequil,  "number of equilibration steps per round");
  argopt_add(ao, "-M",    "%d",   &nstcom,  "interval of centering the structure");
  argopt_add(ao, "-R",    "%b",   &restart, "run a restartable simulation (do not overwrite Zr & rc on output)");
  argopt_add(ao, "-B",    "%b",   &bsim0,   "start a restartable simulation");
  argopt_add(ao, "-C",    "%b",   &norc,    "do not update rc");
  argopt_add(ao, "-W",    "%d",   &nmvtype, "nmove type, 0: Metropolis, 1: heatbath");

  argopt_parse(ao, argc, argv);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  if (norc) updrc = 0;

  nstsave = (int) (nstsav + .5);

  if (rc0 <= 0) rc0 = (real) 0.9;

  /* Metropolis move */
  if (nmvtype == 0) updrc = 0;

  /* make nstcs a multiple of nsted */
  nstcs = (nstcs + nsted - 1) / nsted * nsted;

  if (D >= 15) nstcs = nsted;
  if (nstfb < 0) nstfb = nstcs;

  if (nedxmax < 0) {
    int nedarr[] = {0, 9, 9, 9, 9, 9, 9, 11, 11, 12, 12, 12, 13, 13};
    nedxmax = (D < 14) ? nedarr[D] : 14;
  }

#ifdef MPI
  die_if (nnodes > 1 && neql > 0,
      "no parameter updating with MPI, neql %d\n", neql);
#endif

  if (inode == MASTER) {
    argopt_dump(ao);
  } else { /* change file names for slave nodes */
    if (fninp) fninp = fnappend(fninp, inode);
    if (fnout) fnout = fnappend(fnout, inode);
    if (fnZrtmp) fnZrtmp = fnappend(fnZrtmp, inode);
  }
  argopt_close(ao);
#ifdef MPI
  MPI_Barrier(comm);
#endif
}



typedef struct {
  int nmin, nmax;
  double Z[DG_NMAX + 1]; /* partition function */
  double hist[DG_NMAX + 1]; /* histogram */
  double nup[DG_NMAX + 1][2]; /* number of downward (vertex removing) moves */
  double ndown[DG_NMAX + 1][2]; /* number of upward (vertex adding) moves */
  double Zr[DG_NMAX + 1]; /* Z[n] / Z[n - 1] */
  double Zr1[DG_NMAX + 1]; /* alternative buffer for Zr */
  real rc[DG_NMAX + 1]; /* radius of vertex insertion */
  real rc1[DG_NMAX + 1]; /* alternative buffer of rc */
  double ngpr[DG_NMAX + 1][2]; /* number of good pairs */
  double nedg[DG_NMAX + 1][2]; /* number of edges */
  double ncsp[DG_NMAX + 1][2]; /* no clique separator */
  double fbsm[DG_NMAX + 1][3]; /* hard-sphere weight */
  double ring[DG_NMAX + 1][2]; /* # of ring subgraphs */
  double B[DG_NMAX + 1]; /* virial coefficients */
  double B2[DG_NMAX + 1]; /* virial coefficients, ring route */
  double Bring[DG_NMAX + 1]; /* ring contributions to the virial coefficients */
  double ZZ[DG_NMAX + 1]; /* ratio of Z between two successive pure states */
} gc_t;



/* clear data */
static void gc_cleardata(gc_t *gc)
{
  int i;

  for (i = 0; i <= nmax; i++) {
    gc->hist[i] = 0;
    gc->nup[i][0] = gc->ndown[i][0] = 1e-20;
    gc->nup[i][1] = gc->ndown[i][1] = 0;
    gc->ngpr[i][0] = 1e-20;
    gc->ngpr[i][1] = 0;
    gc->nedg[i][0] = 1e-20;
    gc->nedg[i][1] = 0;
    gc->ncsp[i][0] = 1e-20;
    gc->ncsp[i][1] = 0;
    gc->fbsm[i][0] = 1e-20;
    gc->fbsm[i][1] = gc->fbsm[i][2] = 0;
    gc->ring[i][0] = 1e-20;
    gc->ring[i][1] = 0;
  }
}



INLINE gc_t *gc_open(int nmin, int nmax, real rc0)
{
  gc_t *gc;
  int n;

  xnew(gc, 1);
  gc->nmin = nmin;
  gc->nmax = nmax;

  for (n = 0; n <= nmax; n++) {
    if (nmvtype == 1) { /* heat-bath algorithm */
      if (n == 2) {
        gc->Zr[2] = 0.5;
        gc->rc[2] = 1;
      } else {
        gc->Zr[n] = 1;
        gc->rc[n] = rc0;
      }
    } else { /* Metropolis algorithm */
      if (n <= 2) {
        gc->Zr[n] = 1;
      } else if (n == 3) {
        gc->Zr[3] = Z3rat(D);
      } else if (n == 4) {
        gc->Zr[4] = Z4rat(D) / gc->Zr[3];
      } else {
        gc->Zr[n] = 0.66 * n;
      }
    }
  }
  gc_cleardata(gc);
  if (inode == MASTER)
    loadBring(D, 1, gc->nmax, gc->Bring, fnBring);
#ifdef MPI
  MPI_Bcast(gc->Bring, gc->nmax + 1, MPI_DOUBLE, MASTER, comm);
#endif
  return gc;
}



INLINE void gc_close(gc_t *gc)
{
  free(gc);
}



/* update Zr and possibly rc */
static void gc_update(gc_t *gc, int updrc, double mindata,
    double *Zr, double *rc)
{
  int i, cnt = gc->nmax + 1;
  double r, nZr;

  if (Zr != gc->Zr) memcpy(Zr, gc->Zr, sizeof(Zr[0]) * cnt);
  if (rc != gc->rc) memcpy(rc, gc->rc, sizeof(rc[0]) * cnt);

  for (i = gc->nmin + 1; i <= gc->nmax; i++) {
    if (i == 2) continue; /* exact */
    r = getrrat(gc->nup[i - 1][1], gc->nup[i - 1][0],
        gc->ndown[i][1], gc->ndown[i][0], mindata, 0.9, 1.1);
    if (nmvtype == 1 && updrc && gc->ngpr[i][1] > mindata) {
      /* we treat Zr as the factor used in the acceptance probability
       *  Acc(n -> n-1) = min{1, Zr E}
       * and
       *  Acc(n-1 -> n) = min{1, 1./(Zr E) }
       * So we wish to make Zr <E> ~ 1.
       * The acceptance ratios satisfy
       *  AR(n -> n-1)/AR(n-1 -> n) = rc^D f n (n - 1) / ZZ.
       * where, ZZ = Z(n)/Z(n-1) is the ratio of the actual
       * partition function.
       * Ideally, the ratio of the two AR is 1, so
       *  rc*^D f* n (n - 1) / ZZ = 1. */
      nZr = gc->ngpr[i][0] / gc->ngpr[i][1];
      rc[i] *= pow(r * Zr[i] / nZr, 1. / D);
      Zr[i] = nZr;
    } else { /* update Zr only */
      Zr[i] *= r;
    }
    if (nmvtype == 0) rc[i] = pow(Zr[i], 1. / D);
  }
}



/* compute the partition function
 * assuming gc_update() has been called */
static void gc_computeZ(gc_t *gc, const double *Zr, const double *rc)
{
  int i;
  double z, fbav, fac, nr;

  gc->Z[0] = gc->Z[1] = 1;
  gc->B[0] = gc->B[1] = 1;
  gc->ZZ[0] = gc->ZZ[1] = 1;
  for (fac = 1, z = 1, i = 2; i <= nmax; i++) {
    /* The adjusted acceptance ratios satisfy
     *  AR(n -> n-1) / AR(n-1 -> n) = 1 = rc*^D f* n (n - 1) / ZZ
     * where, ZZ = Z(n)/Z(n-1) is the ratio of the actual
     * partition function. */
    z *= Zr[i];
    if (nmvtype == 1) z *= pow(rc[i], D) * i * (i - 1);
    if (i <= 4) { /* reset the exact values */
      if (i <= 2) z = 1;
      else if (i == 3) {
        printf("%4d: Z3: %.4e to %.4e (exact)\n", inode, z, Z3rat(D));
        z = Z3rat(D);
      } else if (i == 4) {
        printf("%4d: Z4: %.4e vs %.4e (%s)\n", inode, z, Z4rat(D), D > 12 ? "approx." : "exact");
        if (D < 12) z = Z4rat(D);
      }
    }
    fac *= 2./i;
    fbav = gc->fbsm[i][1] / gc->fbsm[i][0];
    gc->Z[i] = z;
    /* B(i)/B(2)^(i-1) = (1 - i) 2^(i-1) / i! Z(i)/Z(2)^(i-1) <fb(i)> */
    gc->B[i] = (1 - i) * fac * z * fbav;
    gc->ZZ[i] = gc->Z[i] / gc->Z[i - 1];
    nr = gc->ring[i][1] / gc->ring[i][0];
    if (nr > 0)
      gc->B2[i] = -fbav / nr * gc->Bring[i];
  }
}



/* print arrays to screen */
static void gc_printZr(const gc_t *gc, const double *Zr, const real *rc)
{
  int i;

  for (i = gc->nmin; i <= gc->nmax; i++)
    printf("%2d %9.6f %8.6f %.6e %12.0f "
        "%.5f %.5f %7.4f %7.4f %10.8f %9.6f %+9.6f %+.6e %+.6e\n",
        i, Zr[i], rc[i], gc->Z[i], gc->hist[i],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0],
        gc->ngpr[i][1] / gc->ngpr[i][0],
        gc->nedg[i][1] / gc->nedg[i][0],
        gc->ncsp[i][1] / gc->ncsp[i][0],
        gc->ring[i][1] / gc->ring[i][0],
        gc->fbsm[i][1] / gc->fbsm[i][0],
        gc->B[i], gc->B2[i]);
}



INLINE double gc_fprintZr(gc_t *gc, FILE *fp,
    const double *Zr, const real *rc)
{
  double tot = 0;
  int i;

  for (i = gc->nmin; i <= gc->nmax; i++) {
    fprintf(fp, "%3d %20.14e %20.14e %18.14f %18.14f "
        "%14.0f %14.0f %14.0f %.14f %.14f ",
        i, gc->ZZ[i], gc->Z[i], Zr[i], rc[i],
        gc->hist[i], gc->nup[i-1][0], gc->ndown[i][0],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0]);
    if (nmvtype == 1)
      fprintf(fp, "%14.0f %.14f ", gc->ngpr[i][0],
          gc->ngpr[i][1] / gc->ngpr[i][0]);
    fprintf(fp, "%14.0f %17.14f ", gc->ring[i][0],
        gc->ring[i][1] / gc->ring[i][0]);
    fprintf(fp, "%14.0f %14.0f %14.0f %18.14f %16.14f %16.14f %+17.14f "
        "%+20.14e %+20.14e\n",
        gc->nedg[i][0], gc->ncsp[i][0], gc->fbsm[i][0],
        gc->nedg[i][1] / gc->nedg[i][0],
        gc->ncsp[i][1] / gc->ncsp[i][0],
        gc->fbsm[i][2] / gc->fbsm[i][0],
        gc->fbsm[i][1] / gc->fbsm[i][0],
        gc->B[i], gc->B2[i]);
    tot += gc->hist[i];
  }
  return tot;
}



#define mkfnZrdef(fn, fndef, d, n, inode) if (fn == NULL) { \
  sprintf(fndef, "Zr%sD%dn%d.dat%d", (nmvtype == 1) ? "h" : "", d, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* save the partition function to file
 * we save the ratio of the successive values, such that one can modify
 * a particular value without affecting later ones */
static int gc_saveZr(gc_t *gc, const char *fn,
    const double *Zr, const real *rc, int dirty)
{
  FILE *fp;
  char fndef[64];
  double tot = 0;

  mkfnZrdef(fn, fndef, D, nmax, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#%s %d %d V4 %d 1 %d %s\n", (nmvtype == 1) ? "H" : "",
      D, gc->nmax, gc->nmin, nedxmax, dirty ? "dirty" : "");
  tot = gc_fprintZr(gc, fp, Zr, rc);
  fclose(fp);
  printf("%4d: saved Zrh to %s, tot %g\n", inode, fn, tot);
  return 0;
}



/* scan data from fp */
static int gc_fscanZr(gc_t *gc, FILE *fp, int loaddata,
    const char *fn, int tag, int ver, double *tot)
{
  int i, n = 0, n0 = gc->nmin, offset, next = 0;
  char s[512], *p;
  double ZZ, Z, rc1;

  *tot = 0;
  offset = (ver == 0) ? 1 : 0;
  if (ver == 0) n0 = 2; else if (ver <= 2) n0 = 3;
  for (i = n0; i <= gc->nmax; i++) {
    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if ( 5 != sscanf(s, "%d%lf%lf%lf%lf%n",
          &n, &ZZ, &Z, &gc->Zr[i], &rc1, &next)
      || i != n + offset) {
      fprintf(stderr, "%s(%d) V%d: ends on line %d(n %d)\n%s", fn, tag, ver, i, n, s);
      break;
    }
    gc->rc[i] = (real) rc1;
    if (nmvtype == 0)
      gc->rc[i] = pow(gc->Zr[i], 1./D);
    if (!loaddata) continue;

    /* try to get additional data */
    p = s + next;
    if ( 5 != sscanf(p, "%lf%lf%lf%lf%lf%n", &gc->hist[i],
              &gc->nup[i-1][0], &gc->ndown[i][0],
              &gc->nup[i-1][1], &gc->ndown[i][1], &next) ) {
      fprintf(stderr, "%s(%d): corrupted(I) on line %d\n%s", fn, tag, i, s);
      break;
    }
    *tot += gc->hist[i];
    if (gc->nup[i-1][0] <= 0) gc->nup[i-1][0] = 1e-20;
    if (gc->ndown[i][0] <= 0) gc->ndown[i][0] = 1e-20;
    gc->nup[i-1][1] *= gc->nup[i-1][0];
    gc->ndown[i][1] *= gc->ndown[i][0];
    p += next;

    if (nmvtype == 1) {
      if ( 2 != sscanf(p, "%lf%lf%n",
                &gc->ngpr[i][0], &gc->ngpr[i][1], &next) ) {
        fprintf(stderr, "%s(%d): corrupted(II) on line %d\n%s", fn, tag, i, s);
        break;
      }
      p += next;
      if (gc->ngpr[i][0] <= 0) gc->ngpr[i][0] = 1e-20;
      gc->ngpr[i][1] *= gc->ngpr[i][0];
    }

    if (ver >= 4) { /* ring data */
      if ( 2 != sscanf(p, "%lf%lf%n",
                &gc->ring[i][0], &gc->ring[i][1], &next) ) {
        fprintf(stderr, "%s(%d): no ring data on line %d\n%s", fn, tag, i, s);
        break;
      }
      p += next;
      if (gc->ring[i][0] <= 0) gc->ring[i][0] = 1e-20;
      gc->ring[i][1] *= gc->ring[i][0];
    }

    if ( 7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf",
           &gc->nedg[i][0], &gc->ncsp[i][0], &gc->fbsm[i][0],
           &gc->nedg[i][1], &gc->ncsp[i][1],
           &gc->fbsm[i][2], &gc->fbsm[i][1]) ) {
      fprintf(stderr, "%s(%d): corrupted(III) on line %d\n%s", fn, tag, i, s);
      break;
    }
    if (gc->nedg[i][0] <= 0) gc->nedg[i][0] = 1e-20;
    if (gc->ncsp[i][0] <= 0) gc->ncsp[i][0] = 1e-20;
    if (gc->fbsm[i][0] <= 0) gc->fbsm[i][0] = 1e-20;
    gc->nedg[i][1] *= gc->nedg[i][0];
    gc->ncsp[i][1] *= gc->ncsp[i][0];
    gc->fbsm[i][1] *= gc->fbsm[i][0];
    gc->fbsm[i][2] *= gc->fbsm[i][0];
  }
  /* should be gc->nmax + 1 if successful */
  return i;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int gc_loadZr(gc_t *gc, const char *fn, int loaddata)
{
  FILE *fp;
  int i = -1, d = 0, n = 0, n0 = 1, next = 0, ver = 0;
  char s[512], fndef[64], sver[16] = "";
  double tot;

  mkfnZrdef(fn, fndef, D, nmax, inode);
  xfopen(fp, fn, "r", return -1);

  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if ( s[0] != '#' || (nmvtype == 1 && s[1] != 'H')
    || sscanf(s + 2, "%d%d%8s%d%n", &d, &n, sver, &n0, &next) != 4
    || d != D) {
    fprintf(stderr, "%s: D %d vs %d\n%s", fn, d, D, s);
    fclose(fp);
    return -1;
  }
  ver = atoi(sver[0] == 'V' ? sver + 1 : sver); /* strip V */
  if (loaddata && (ver <= 2 || strstr(s, "dirty") || strstr(s, "DIRTY")) ) {
    fprintf(stderr, "%s: %s has no restartable data\n", fn, sver);
    loaddata = 0;
  }
  i = gc_fscanZr(gc, fp, loaddata, fn, -1, ver, &tot);
  fclose(fp);
  printf("%4d: loaded Zr %s, V%d, tot %g, data %d\n",
      inode, fn, ver, tot, loaddata);
  return i;
}



/* accumulate data */
static void gc_accumdata(gc_t *gc, const dg_t *g, double t,
    int nstcs, int nstfb)
{
  int ned = -1, n = g->n, ncs, err, errnr;
  double fb = 0, sc;
  double nr = 0;
  static int degs[DG_NMAX], nlookup = 0;

  if (nlookup <= 0)
    /* we will always compute fb if a lookup table is available
     * i.e., for a small n. But in the n = 8 case, initializing
     * the lookup table would cost 30s-70s,
     * which is a bit long for a short simulation */
    nlookup = (nsteps < 1e9) ? 7 : DGMAP_NMAX;

  ned = dg_degs(g, degs);
  gc->nedg[n][0] += 1;
  gc->nedg[n][1] += ned;

  if (n <= nlookup || ((int) fmod(t + .5, nstcs) == 0) ) {
    /* check if the graph has a clique separator
     * but in special cases, fb and nr are computed as well */
    if (n <= 3) { /* assuming biconnectivity */
      static double fbarr[4] = {1, 1, -1, -1};
      err = errnr = 0;
      ncs = 1;
      fb = fbarr[n];
      nr = 1;
    } else if (n <= nlookup) { /* lookup table */
      code_t code;
      unqid_t uid;
      dgmap_t *m = dgmap_ + n;

      dg_encode(g, &code);
      ncs = (dg_ncsep_lookuplow(g, code) == 0);
      dgmap_init(m, n);
      uid = m->map[code];
      fb = dg_hsfb_lookuplow(n, uid);
      nr = dg_nring_lookuplow(n, uid);
      err = errnr = 0;
    } else {
      /* this function implicitly computes the clique separator
       * with very small overhead */
      sc = dg_rhsc_spec0(g, 0, &ned, degs, &err);
      if (err == 0) {
        ncs = (fabs(sc) > 1e-3);
        fb = DG_SC2FB(sc, ned);
      } else { /* no clique separator */
        ncs = 1;
      }
      errnr = 1; /* indicate we don't have nr */
    }
    gc->ncsp[n][0] += 1;
    gc->ncsp[n][1] += ncs;

    if (n <= nlookup || (nstfb > 0 && (int) fmod(t + .5, nstfb) == 0)) {
      /* compute fb, if it is cheap */
      if ( err ) { /* if dg_rhsc_spec0() fails, no clique separator */
        if (ned > n + nedxmax && n > nedxmax + 2) {
          fprintf(stderr, "%d: t %.6e/%g D %d assume fb = 0 n %d ned %d, del %d\n",
              inode, t, nsteps, D, n, ned, ned - n);
          fb = 0;
        } else {
          fb = dg_hsfb_mixed0(g, 1, &ned, degs);
        }
      }
      gc->fbsm[n][0] += 1;
      gc->fbsm[n][1] += fb;
      gc->fbsm[n][2] += fabs(fb);

      /* compute the ring content */
      if (errnr) {
        nr = dg_nring_spec0(g, &ned, degs, &errnr);
        if (errnr) nr = dg_nring_direct(g);
      }
      gc->ring[n][0] += 1;
      gc->ring[n][1] += nr;
    }
  }
}



/* get a pair of vertices within rc */
INLINE int getpair(int *pi, int *pj, const dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int i, j, npr = 0, id, n = g->n;
  static int pr[DG_NMAX*DG_NMAX];
  code_t mask = mkbitsmask(n);

  for (i = 0; i < n; i++) { /* the vertex to remove */
    if ( n > 2 && !dg_biconnectedvs(g, mask ^ MKBIT(i)) ) continue;
    for (j = 0; j < n; j++) { /* root */
      if ( j == i ) continue;
      if ( r2ij[i][j] < rc * rc )
        pr[npr++] = i * DG_NMAX + j;
    }
  }

  if (npr <= 0) return 0;
  if (pi != NULL && pj != NULL) {
    id = pr[ (int) (rnd0() * npr) ];
    *pi = id / DG_NMAX;
    *pj = id - (*pi) * DG_NMAX;
  }
  return npr;
}



/* grand canonical simulation, main function */
static int mcgc(int nmin, int nmax, double nsteps, double mcamp,
    int neql, double nequil, int nstsave)
{
  double t, cacc = 0, ctot = 0;
  int ninit, i, j, k, n, npr, deg, acc, it, ieql, loaddata;
  gc_t *gc;
  dg_t *g, *ng;
  rvn_t x[DG_NMAX], xi;
  real r2ij[DG_NMAX][DG_NMAX], r2i[DG_NMAX];

  gc = gc_open(nmin, nmax, rc0);
  loaddata = restart && !bsim0;
  gc_loadZr(gc, fninp, loaddata);
#ifdef MPI
  /* make sure the same parameters */
  MPI_Bcast(gc->Zr, gc->nmax + 1, MPI_DOUBLE, MASTER, comm);
  MPI_Bcast(gc->rc, gc->nmax + 1, MPI_MYREAL, MASTER, comm);
  if (nnodes > 1) /* scramble the random number generator */
    mtscramble(inode * 2034091783u + time(NULL));
  MPI_Barrier(comm);
  if (inode == MASTER)
    for (i = nmin; i <= nmax; i++)
      printf("%3d: %.6f %.6f\n", i, gc->Zr[i], gc->rc[i]);
  MPI_Barrier(comm);
#else
  for (i = nmin; i <= nmax; i++)
    printf("%3d: %.6f %.6f\n", i, gc->Zr[i], gc->rc[i]);
#endif

  /* we limit the initial # of vertices to 8 for larger
   * fully-connected diagrams are harder to equilibrate */
  if (nmin < 8) {
    int nn = (8 > nmax) ? nmax : 8;
    ninit = nmin + (int) (rnd0() * (nn - nmin));
  } else
    ninit = nmin;
  g = dg_open(nmax);
  ng = dg_open(nmax);
  initx(x, nmax);
  /* initial n in nmin */
  calcr2ij(r2ij, x, ninit);
  mkgraphr2ij(g, r2ij, 1, ninit);
  printf("%4d/%-4d: n %d [%d, %d]\n", inode, nnodes, g->n, nmin, nmax);

  ieql = (neql > 0) ? 1 : neql; /* stage of equilibrations if > 0 */
  /* ieql < 0 means no data collection */

  for (it = 1, t = 1; ieql || t <= nsteps; t += 1, it++) {
    n = g->n;
    die_if (n < nmin || n > nmax, "bad n %d, t %g\n", n, t);
    if ( rnd0() < ratn ) { /* switching the ensemble, O(n) */
      if ( rnd0() < 0.5 ) { /* add an vertex */
        if ( n >= nmax ) goto STEP_END;
        /* attach a new vertex to a random vertex */
        i = (n == 1) ? 0 : (int) (rnd0() * n);
        rvn_inc(rvn_rndball(x[n], gc->rc[n + 1]), x[i]);
        /* test if adding x[n] leaves the diagram biconnected */
        for (deg = 0, j = 0; j < n; j++)
          if ( (r2i[j] = rvn_dist2(x[n], x[j])) < 1 )
            ++deg;
        if ( n == 1 ) deg *= 2;
        /* the extended configuration is biconnected
         * if x[n] is connected two vertices */
        if ( ieql >= 0 ) gc->nup[n][0] += 1;
        if ( deg >= 2 ) {
          g->n = n + 1;
          for (j = 0; j < n; j++) {
            r2ij[j][n] = r2ij[n][j] = r2i[j];
            if ( r2i[j] < 1)
              dg_link(g, j, n);
            else
              dg_unlink(g, j, n);
          }
          if ( nmvtype == 1 ) { /* heat-bath algorithm */
            npr = getpair(NULL, NULL, g, r2ij, gc->rc[n + 1]);
            if (rnd0() < 1./(gc->Zr[n + 1] * npr)) {
              /* accept the move */
              if ( ieql >= 0 ) gc->nup[n][1] += 1;
            } else { /* recover */
              dg_remove1(g, g, n);
            }
          } else { /* Metropolis algorithm */
            if ( ieql >= 0 ) gc->nup[n][1] += 1;
          }
        }
      } else if ( n > nmin ) { /* remove a random vertex */
        if ( ieql >= 0 ) gc->ndown[n][0] += 1;
        if ( nmvtype == 1 ) { /* heat-bath removal */
          if ( (npr = getpair(&i, &j, g, r2ij, gc->rc[n])) <= 0 )
            goto STEP_END;
          if ( ieql >= 0 ) {
            gc->ngpr[n][0] += 1;
            gc->ngpr[n][1] += npr;
          }
          die_if (n == 2 && npr != 2, "n %d, npr %d\n", n, npr);
          if (rnd0() >= npr * gc->Zr[n]) goto STEP_END;
        } else { /* Metropolis removal */
          i = randpair(n, &j);
          die_if (i == j || i >= n || j >= n, "bad i %d, j %d, n %d\n", i, j, g->n);
          /* test if * the graph is biconnected without i
           * and the pair i and j are connected */
          if ( ( n > 2
              && !dg_biconnectedvs(g, mkbitsmask(n) ^ MKBIT(i)) )
           || rvn_dist2(x[i], x[j]) >= dblsqr(gc->rc[n]) )
            goto STEP_END;
        }

        /* accepted the removal */
        if ( ieql >= 0 ) gc->ndown[n][1] += 1;
        dg_remove1(g, g, i); /* n is decreased by 1 here */
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
      }
    }
    else /* configuration sampling */
    {
      if (n == 1) acc = 1;
      else BCSTEPR2(acc, i, n, g, ng, x, xi, r2ij, r2i, mcamp, gaussdisp);
      cacc += acc;
      ctot += 1;
    }
STEP_END:
    if ( ieql >= 0 ) {
      gc->hist[g->n] += 1;
      if (it % nsted == 0)
        gc_accumdata(gc, g, t, nstcs, nstfb);
    }
    if (it % nstcom == 0) /* remove the center of mass motion */
      shiftr2ij(g, x, r2ij);
    if (ieql != 0) { /* equilibration */
      if ((int) fmod(t + .5, nequil) == 0) {
        if (ieql > 0) { /* update parameters */
          gc_update(gc, updrc, mindata, gc->Zr, gc->rc);
          gc_computeZ(gc, gc->Zr, gc->rc);
          gc_saveZr(gc, fnZrtmp, gc->Zr, gc->rc, 1);
          gc_printZr(gc, gc->Zr, gc->rc);
          printf("D %d, equilibration stage %d/%d, t %g\n", D, ieql, neql, t);
          gc_cleardata(gc);
          if (++ieql > neql) /* stop equilibration */
            ieql = 0;
        } else {
          printf("%d, D %d, t %g: equilibrated, n %d\n", inode, D, t, g->n);
          ieql = 0;
        }
        t = 0; it = 0; /* reset time */
      }
    } else { /* production */
      if (it % nstsave == 0 || t > nsteps - 0.5) {
        /* compute node-specific Z */
        gc_update(gc, updrc, mindata, gc->Zr1, gc->rc1);
        gc_computeZ(gc, gc->Zr1, gc->rc1);
        if (restart) {
          gc_saveZr(gc, fnout, gc->Zr, gc->rc, 0);
          if (inode == MASTER)
            gc_saveZr(gc, fnZrtmp, gc->Zr1, gc->rc1, 1);
        } else {
          gc_saveZr(gc, fnout, gc->Zr1, gc->rc1, 1);
        }
        it = 0;
        if (inode == MASTER) mtsave(NULL);
      }
    }
  }

  printf("%d: cacc %g\n", inode, cacc/ctot);
  if (inode == MASTER)
    gc_printZr(gc, gc->Zr1, gc->rc1);
  dg_close(g);
  dg_close(ng);
  return 0;
}



int main(int argc, char **argv)
{
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnodes);
#endif
  doargs(argc, argv);
  mcgc(nmin, nmax, nsteps, mcamp, neql, nequil, nstsave);
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

