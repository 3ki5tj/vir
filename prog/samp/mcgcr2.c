/* grand-canonical simulation to compute the partition function
 * improved version for higher dimensions */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dgrjw.h"
#include "dgring.h"
#include "mcutil.h"



int nmin = 1; /* the minimal order of virial coefficients */
int nmax = DG_NMAX; /* the maximal order of virial coefficients */
int mtiers = 2;
real mcamp = 1.5f;
double nequil = 10000000;
double nsteps = 10000000;
double ratn = 0.5;
int gaussdisp = 0;
real rc0 = -1;
real sr0 = 1;
double mindata = 100; /* minimal # of data points to make update */
int updrc = 0;

/* common and shared file names */
char *fninp0, *fnZrr0, *fnZr0, *fnZrrtmp = "Zrr.tmp";

/* thread-dependent file names */
char *fninp, *fnZrr, *fnZr;
#pragma omp threadprivate(fninp, fnZrr, fnZr)

char *fnBring = NULL; /* ring contribution file */

unsigned rngseed = 0;

int nstcom = 10000;
int nsted = 10;
int nstcs = -1;
int nstfb = -1;
/* maximal number of extra edges to invoke fb calculation
 * reference T60 timing data (no clique separator):
 * nedxmax  time per evaluation (ms)
 *  10         15
 *  11         50
 *  12        150
 * */
int nedxmax = -1;

int neql = -1; /* number of equilibration stages
                 == 0: disable
                 >  0: Zr, sr are updated after each round
                 <  0: Zr, sr are not updated, one round */
int nstsave = 1000000000; /* interval of saving data */
int nstcheck = 0;

int restart = 0; /* restartable simulation */
int bsim0 = 0; /* first simulation */

/* nmove type, 0: Metropolis, 1: heat-bath */
int nmvtype = 1;



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int m1 = 1;
  double nstsav = nstsave;

  ao = argopt_open(0);
  argopt_add(ao, "-z",    "%d",   &nmin,    "minimal order of virial coefficients");
  argopt_add(ao, "-n",    "%d",   &nmax,    "maximal order of virial coefficients");
  argopt_add(ao, "-m",    "%d",   &m1,      "number of intermediate ensembles per order of virial coefficients");
  argopt_add(ao, "-1",    "%lf",  &nsteps,  "number of simulation steps");
  argopt_add(ao, "-a",    "%r",   &mcamp,   "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-r",    "%lf",  &ratn,    "rate of particle moves");
  argopt_add(ao, "-c",    "%r",   &rc0,     "radius of particle insersion");
  argopt_add(ao, "-s",    "%r",   &sr0,     "default scaling of distance between the constrained vertex pair");
  argopt_add(ao, "-i",    NULL,   &fninp0,  "input file");
  argopt_add(ao, "-o",    NULL,   &fnZrr0,  "output file");
  argopt_add(ao, "-O",    NULL,   &fnZr0,   "output file (compact)");
  argopt_add(ao, "-T",    NULL,   &fnZrrtmp,"temporary output file");
  argopt_add(ao, "-I",    NULL,   &fnBring, "ring contribution file");
  argopt_add(ao, "-G",    "%d",   &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-A",    "%d",   &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-F",    "%d",   &nstfb,   "interval of computing fb");
  argopt_add(ao, "-X",    "%d",   &nedxmax, "maximal number of extra edges for computing fb");
  argopt_add(ao, "-Q",    "%lf",  &mindata, "minimal number of data points to warrant parameter update");
  argopt_add(ao, "-C",    "%b",   &updrc,   "update rc");
  argopt_add(ao, "-q",    "%lf",  &nstsav,  "interval of saving data");
  argopt_add(ao, "-E",    "%d",   &neql,    "number of equilibration stages, 0: disable, > 0: update Zr & sr, < 0: 1 round no update");
  argopt_add(ao, "-0",    "%lf",  &nequil,  "number of equilibration steps per round");
  argopt_add(ao, "-M",    "%d",   &nstcom,  "interval of centering the structure");
  argopt_add(ao, "-K",    "%d",   &nstcheck,"interval of checking");
  argopt_add(ao, "-R",    "%b",   &restart, "run a restartable simulation");
  argopt_add(ao, "-B",    "%b",   &bsim0,   "start a restartable simulation");
  argopt_add(ao, "-W",    "%d",   &nmvtype, "nmove type, 0: Metropolis, 1: heatbath");
  argopt_add(ao, "--rng", "%u",   &rngseed, "set the RNG seed offset");

  argopt_parse(ao, argc, argv);

  die_if (nmax > DG_NMAX, "too many %d atoms >= %d\n", nmax, DG_NMAX);

  mtiers = m1 + 1;
  die_if (mtiers <= 1, "mtiers %d should be at least 2\n", mtiers);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  nstsave = (int) (nstsav + .5);

  if (rc0 <= 0) {
    if (nmvtype == 1) {
      /* this value works for large n */
      rc0 = (real) (1. / mtiers);
    } else {
      rc0 = (real) (4. / sqrt(D));
      if (rc0 > 1) rc0 = 1;
    }
  }

  if (nstcs >= 0) /* round nstcs to a multiple of nsted */
    nstcs = (nstcs + nsted - 1) / nsted * nsted;
  if (D >= 15 && nstcs < 0) nstcs = nsted;
  if (nstfb < 0) nstfb = nstcs;

  if (nedxmax < 0) {
    int nedarr[] = {0, 9, 9, 9, 9, 9, 9, 11, 11, 12, 12, 12, 13, 13};
    nedxmax = (D < 14) ? nedarr[D] : 14;
  }

  die_if (nnodes > 1 && neql > 0,
      "no parameter updating with MPI, neql %d\n", neql);

  if (inode == MASTER) {
    argopt_dump(ao);
  }
  argopt_close(ao);
#ifdef MPI
  MPI_Barrier(comm);
#endif
}



/* append node ids after the file names */
INLINE void appendfns(const char *fninp0, const char *fnZrr0,
    const char *fnZr0)
{
  if (inode != MASTER) {
    if (fninp0) fninp = fnappend(fninp0,    inode);
    if (fnZrr0) fnZrr = fnappend(fnZrr0,    inode);
    if (fnZr0)  fnZr  = fnappend(fnZr0,     inode);
  } else {
    if (fninp0) fninp = ssdup(fninp0);
    if (fnZrr0) fnZrr = ssdup(fnZrr0);
    if (fnZr0)  fnZr  = ssdup(fnZr0);
  }
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
  double (*nup)[2]; /* number of downward (vertex removing) moves */
  double (*ndown)[2]; /* number of upward (vertex adding) moves */
  double *Zr; /* Z[iens] / Z[iens - 1] */
  real *rc; /* radius of the vertex insertion */
  real *sr; /* scaling between the special pair of vertices */
  double *Zr1; /* alternative buffer for Zr */
  real *rc1, *sr1; /* alternative buffer for rc and sr */
  double (*nedg)[2]; /* number of edges */
  double (*ncsp)[2]; /* no clique separator */
  double (*fbsm)[3]; /* hard-sphere weight */
  double (*ring)[2]; /* # of ring subgraphs */
  double *B; /* virial coefficients */
  double *B2; /* virial coefficients, ring route */
  double *Bring; /* ring contributions to the virial coefficients */
  double *Zn; /* partition function at pure states */
  double *ZZ; /* ratio of Z between two successive pure states */
} gc_t;

#ifndef MPI
/* shared variables for OpenMP */
double *gc_Zr;
real *gc_rc;
real *gc_sr;
double *gc_Bring;
#endif


/* clear data */
static void gc_cleardata(gc_t *gc)
{
  int i, n;

  for (i = 0; i < gc->nens; i++) {
    gc->hist[i] = 0;
    gc->nup[i][0] = gc->ndown[i][0] = 1e-20;
    gc->nup[i][1] = gc->ndown[i][1] = 0;
  }
  for (n = 0; n <= gc->nmax; n++) {
    gc->nedg[n][0] = 1e-20;
    gc->nedg[n][1] = 0;
    gc->ncsp[n][0] = 1e-20;
    gc->ncsp[n][1] = 0;
    gc->fbsm[n][0] = 1e-20;
    gc->fbsm[n][1] = gc->fbsm[n][2] = 0;
    gc->ring[n][0] = 1e-20;
    gc->ring[n][1] = 0;
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
  xnew(gc->nedg, nmax + 1);
  xnew(gc->ncsp, nmax + 1);
  xnew(gc->fbsm, nmax + 1);
  xnew(gc->ring, nmax + 1);
  xnew(gc->B, nmax + 1);
  xnew(gc->B2, nmax + 1);
  xnew(gc->Bring, nmax + 1);
  xnew(gc->Zn, nmax + 1);
  xnew(gc->ZZ, nmax + 1);

  /* temporary setting */
  if (m > 2 && rc0 >= 1) rc0 = 0.8;
  for (i = 0; i < nensmax; i++) {
    gc->n[i] = (i + m - 1) / m;
    /* type == 0 mean GCX_PURE */
    gc->type[i] = (i % m);
  }

  for (i = 0; i < nensmax; i++) {
    n = gc->n[i];

    if ( gc->type[i] == GCX_PURE ) {
      rc = 1;
    } else {
      /* rc = pow(rc0, 1.*(m - i % m)/(m - 1)); */
      /* linear interpolation, including the boundary value 1
       * By default, rc0 = 1/m. When i == 1, rc = rc0 = 1/m
       * when i == m - 1, rc = (m - 1)/m */
      rc = rc0 + (1. - rc0) * (i % m - 1) / (m - 1);
    }
    gc->rc[i] = (real) rc;

    if (gc->type[i] == GCX_PURE && n >= 1)
      gc->Zr[i] = 1.0/n; /* inverse number of nonarticulated pairs */
    else gc->Zr[i] = 1;
    gc->Z[i] = 1; /* absolute partition function */
  }

  /* sr[i] requires a separate loop for it needs rc[i + 1] */
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

  for (n = 0; n <= nmax; n++)
    gc->ZZ[n] = gc->Zn[n] = gc->B[n] = gc->B2[n] = 1;
  gc_cleardata(gc);
  if (inode == MASTER)
    loadBring(D, 1, gc->nmax, gc->Bring, fnBring);
#ifdef MPI
  MPI_Bcast(gc->Bring, gc->nmax + 1, MPI_DOUBLE, MASTER, comm);
#else
  if (inode == MASTER) {
    gc_Zr = gc->Zr;
    gc_rc = gc->rc;
    gc_sr = gc->sr;
    gc_Bring = gc->Bring;
  }
  /* wait the master */
  #pragma omp barrier
  if (inode != MASTER)
    for (n = 0; n <= gc->nmax; n++)
      gc->Bring[n] = gc_Bring[n];
#endif
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
  free(gc->nedg);
  free(gc->ncsp);
  free(gc->fbsm);
  free(gc->ring);
  free(gc->B);
  free(gc->B2);
  free(gc->Bring);
  free(gc->Zn);
  free(gc->ZZ);
  free(gc);
}



/* accumulate data */
static void gc_accumdata(gc_t *gc, const dg_t *g, double t,
    int nstcs, int nstfb)
{
  int ned = -1, n = g->n, ncs, err, errnr;
  double sc, fb = 0;
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

  if (nstcs <= 0) return;
  if (n <= nlookup || rnd0() * nstcs / nsted < 1) {
    /* check if the graph has a clique separator
     * but in special cases, fb and nr are computed as well */
    if (n <= 3) { /* assuming biconnectivity */
      static double fbarr[4] = {1, 1, -1, -1};
      err = errnr = 0;
      ncs = 1;
      fb = fbarr[n];
      nr = 1;
    } else if (n <= nlookup) { /* lookup table */
      dgword_t code;
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
       * with a very small overhead */
      sc = dg_rhsc_spec0(g, 0, 1, &ned, degs, &err);
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

    if (nstfb <= 0) return;
    if (n <= nlookup || rnd0() * nstfb/nstcs < 1) {
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



/* update parameters */
static void gc_update(gc_t *gc, double mindata, int updrc,
    double *Zr, real *rc, real *sr)
{
  int i;
  double r;

  if (Zr != gc->Zr) memcpy(Zr, gc->Zr, sizeof(Zr[0]) * gc->nens);
  if (rc != gc->rc) memcpy(rc, gc->rc, sizeof(rc[0]) * gc->nens);
  if (sr != gc->sr) memcpy(sr, gc->sr, sizeof(sr[0]) * gc->nens);
  for (i = gc->ens0 + 1; i < gc->nens; i++) {
    r = getrrat(gc->nup[i-1][1], gc->nup[i-1][0],
        gc->ndown[i][1], gc->ndown[i][0], mindata, 0.9, 1.1);
    if (gc->type[i - 1] == GCX_PURE) {
      Zr[i] *= r;
    } else {
      /* A_down / A_up = s^D / Zr (1)
       * ideally 1 = A*_down / A*_up = s*^D / Zr (2)
       * (2)/(1) --> A_up / A_down = (s* / s)^D */
      sr[i - 1] *= pow(r, 1./D);
    }
  }

  if ( !updrc ) return;

  /* update rc, use a separate loop to be clean */
  for (i = gc->ens0; i < gc->nens - 1; i++) {
    double r1, r2, r3, r4, newrc;

    /* pure states don't have rc */
    if ( gc->type[i] == GCX_PURE || gc->n[i] <= 2 ) continue;

    /* since updating rc is dangerous, we do it cautiously */
    if ( gc->nup[i-1][1] < mindata || gc->ndown[i][1] < mindata
      || gc->nup[i][1] < mindata || gc->ndown[i+1][1] < mindata )
      continue;
    r1 = gc->nup[i-1][1] / gc->nup[i-1][0];
    r2 = gc->ndown[i][1] / gc->ndown[i][0];
    r3 = gc->nup[i][1] / gc->nup[i][0];
    r4 = gc->ndown[i+1][1] / gc->ndown[i+1][0];

/* the following parameters are defined as macros
 * so that they can be modified on compiling time */
#ifndef GCFLATNESS
#define GCFLATNESS 0.05
#endif
#ifndef GCEXPONENT
/* maybe needs to be larger for larger D */
#define GCEXPONENT (4.0/D)
#endif
#ifndef GCROUND
/* this number needs to be smaller for larger D */
#define GCROUND 0.001
#endif

    /* make sure the transition rates between i-1 and i are equal
     * those between i and i+1 are also equal */
    if ( fabs(r1 - r2) > GCFLATNESS * (r1 + r2)
      || fabs(r3 - r4) > GCFLATNESS * (r3 + r4) )
      continue;
    r1 = (r1 + r2)/2;
    r3 = (r3 + r4)/2;
    /* multiply rc by a moderate power of r1/r3
     * limit the amount of updating within 10% */
    r = dblconfine(pow(r1 / r3, GCEXPONENT), 0.9, 1.1);
    /* we round rc to a multiple of 0.01 */
    newrc = dblconfine(dblround(r * rc[i], GCROUND), 0.05, 1.0);
    r = newrc / rc[i];
    if (fabs(r - 1) > 1e-6) {
      printf("%3d: iens %3d, n %2d, updating rc %.3f -> %.3f\n",
        inode, i, gc->n[i], rc[i], newrc);
      rc[i] = newrc;
      sr[i] /= r;
      if (gc->type[i - 1] != GCX_PURE)
        sr[i - 1] *= r;
    }
  }
}



/* compute the partition function and virial coefficients */
static void gc_computeZ(gc_t *gc,
    const double *Zr, const real *rc, const real *sr)
{
  double fac = 1, z = 1, fbav, nr, prevZ = 1;
  int n, i;

  if (Zr == NULL) Zr = gc->Zr;
  if (rc == NULL) rc = gc->rc;
  if (sr == NULL) sr = gc->sr;
  for (i = 1; i < gc->nens; i++) {
    n = gc->n[i];
    /* compute the partition function */
    /* Note: Zr[i] == 1 unless i-1 or i is a pure state */
    if (gc->type[i - 1] == GCX_PURE) {
      z *= Zr[i] * pow(rc[i], D);
    } else {
      /* Zr[i] != 1 only if i is a pure state */
      z *= Zr[i] * pow(sr[i - 1], D);
    }
    /* use the the exact result if available */
    if (gc->type[i] == GCX_PURE) {
      if (nmvtype == 1)
        z *= n * (n - 1) * .5; /* pair generation bias */
      if (n == 3 && nmin < 3)
        printf("%4d: Z3: %.4e vs %.4e (exact)\n", inode, z, Z3rat(D));
      if (n == 4)
        printf("%4d: Z4: %.4e vs %.4e (%s)\n", inode, z, Z4rat(D), (D > 12) ? "approx." : "exact");
      if (n <= 2) z = 1;
      else if (n == 3) z = Z3rat(D);
      else if (n == 4 && D <= 12) z = Z4rat(D);
    }
    gc->Z[i] = z;
    if (gc->type[i] == GCX_PURE) {
      if (n >= 2) fac *= 2./n;
      fbav = (n <= 2) ? 1 : gc->fbsm[n][1] / gc->fbsm[n][0];
      /* Bn/B2^(n-1) = (1 - n) 2^(n - 1) / n! Zn/Z2^(n-1) < fb > */
      gc->B[n] = (1. - n) * fac * z * fbav;
      gc->Zn[n] = gc->Z[i];
      gc->ZZ[n] = gc->Zn[n] / prevZ;
      prevZ = gc->Zn[n];
      nr = gc->ring[n][1] / gc->ring[n][0];
      if (nr > 0)
        gc->B2[n] = -fbav / nr * gc->Bring[n];
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
        printf("%6.3f %8.6f %+9.7f %+9.7f %+.6e %+.6e",
          gc->nedg[n][1] / gc->nedg[n][0],
          gc->ncsp[n][1] / gc->ncsp[n][0],
          gc->ring[n][1] / gc->ring[n][0],
          gc->fbsm[n][1] / gc->fbsm[n][0],
          gc->B[n], gc->B2[n]);
      }
    }
    printf("\n");
  }
}



#define mkfnZrdef(fn, fndef, d, m, n, inode) if (fn == NULL) { \
  sprintf(fndef, "ZrD%dr%dn%d.dat%d", d, m, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* save a compact (pure-state) Zr file, mimic the output of mcgc.c
 * this file is also gnuplot-friendly */
static int gc_saveZr(gc_t *gc, const char *fn)
{
  FILE *fp;
  int i, n;
  char fndef[64];
  double tot = 0;

  mkfnZrdef(fn, fndef, D, gc->m - 1, gc->nmax, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d V4 %d 1 %d\n", D, gc->nmax,
      gc->nmin, nedxmax);
  for (i = gc->ens0; i < gc->nens; i++) {
    if (gc->type[i] != GCX_PURE) continue;
    n = gc->n[i];
    fprintf(fp, "%3d %18.14f %20.14e %14.0f %14.0f %14.0f "
        "%.14f %.14f %14.0f %17.14f %14.0f %14.0f %14.0f "
        "%18.14f %16.14f %16.14f %+17.14f %+20.14e %+20.14e\n",
        n, gc->ZZ[n], gc->Z[i], gc->hist[i],
        gc->nup[i-1][0], gc->ndown[i][0],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0],
        gc->ring[n][0], gc->ring[n][1] / gc->ring[n][0],
        gc->nedg[n][0], gc->ncsp[n][0], gc->fbsm[n][0],
        gc->nedg[n][1] / gc->nedg[n][0], gc->ncsp[n][1] / gc->ncsp[n][0],
        gc->fbsm[n][2] / gc->fbsm[n][0], gc->fbsm[n][1] / gc->fbsm[n][0],
        gc->B[n], gc->B2[n]);
    tot += gc->hist[i];
  }
  fclose(fp);
  printf("%4d: saved compact Zr to %s, total %g\n", inode, fn, tot);
  return 0;
}



/* print Zrr data to file */
INLINE double gc_fprintZrr(gc_t *gc, FILE *fp,
    const double *Zr, const real *rc, const real *sr)
{
  double tot = 0;
  int i;

  for (i = 1; i < gc->nens; i++) {
    fprintf(fp, "%4d %d %3d %20.14e %20.14e %18.14f %18.14f "
        "%14.0f %14.0f %14.0f %.14f %.14f ",
        i, gc->type[i], gc->n[i], Zr[i], gc->Z[i], rc[i], sr[i],
        gc->hist[i], gc->nup[i-1][0], gc->ndown[i][0],
        gc->nup[i-1][1] / gc->nup[i-1][0],
        gc->ndown[i][1] / gc->ndown[i][0]);
    if (gc->type[i] == GCX_PURE) {
      int n = gc->n[i];
      fprintf(fp, "%14.0f %17.14f ", gc->ring[n][0],
        gc->ring[n][1] / gc->ring[n][0]);
      fprintf(fp, "%14.0f %14.0f %14.0f ",
        gc->nedg[n][0], gc->ncsp[n][0], gc->fbsm[n][0]);
      fprintf(fp, "%18.14f %16.14f %16.14f %+17.14f %+20.14e %+20.14e ",
        gc->nedg[n][1] / gc->nedg[n][0],
        gc->ncsp[n][1] / gc->ncsp[n][0],
        gc->fbsm[n][2] / gc->fbsm[n][0], /* < |fb| > */
        gc->fbsm[n][1] / gc->fbsm[n][0], /* <  fb  > */
        gc->B[n], gc->B2[n]);
    }
    fprintf(fp, "\n");
    tot += gc->hist[i];
  }
  return tot;
}



#define mkfnZrrdef(fn, fndef, d, m, n, inode) if (fn == NULL) { \
  sprintf(fndef, "ZrrD%dr%dn%d.dat%d", d, m, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* save all data to file */
static int gc_saveZrr(gc_t *gc, const char *fn,
    const double *Zr, const real *rc, const real *sr, int dirty)
{
  FILE *fp;
  char fndef[64];
  double tot = 0;

  mkfnZrrdef(fn, fndef, D, gc->m - 1, gc->nmax, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#%c %d %d %d %d %d %d V4 1 %d %s\n", nmvtype ? 'R' : 'S',
      D, gc->nens, gc->ens0, gc->nmin, gc->nmax, gc->m, nedxmax,
      dirty ? "dirty" : "");
  tot = gc_fprintZrr(gc, fp, Zr, rc, sr);
  fclose(fp);
  printf("%4d: saved Zrr to %s, tot %g\n", inode, fn, tot);
  return 0;
}



/* scan data from fp */
static int gc_fscanZrr(gc_t *gc, FILE *fp, int loaddata,
    const char *fn, int tag, int ver, double *tot)
{
  int i, d = 0, tp = 0, n = 0, next = 0;
  char s[512], *p;
  double Zr, Z, rc, sr;

  *tot = 0;
  for (i = 1; i < gc->nens; i++) {
    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if (7 != sscanf(s, "%d%d%d%lf%lf%lf%lf%n",
                       &d, &tp, &n, &Zr, &Z, &rc, &sr, &next)
        || i != d || tp != gc->type[i] || n != gc->n[i]) {
      fprintf(stderr, "%s(%d) ends on line %d\n%s", fn, tag, i, s);
      break;
    }
    if (Zr >= 0) gc->Zr[i] = Zr;
    if (rc >= 0) gc->rc[i] = rc;
    if (sr >= 0) {
      if (ver >= 2) gc->sr[i] = sr;
      else gc->sr[i - 1] = sr; /* old version */
    }

    if (!loaddata) continue;
    /* try to get additional data */
    p = s + next;
    if (5 != sscanf(p, "%lf%lf%lf%lf%lf%n",
          &gc->hist[i], &gc->nup[i-1][0], &gc->ndown[i][0],
          &gc->nup[i-1][1], &gc->ndown[i][1], &next) ) {
      fprintf(stderr, "%s(%d) corrupted on line %d\n%s", fn, tag, i, s);
      break;
    }
    *tot += gc->hist[i];
    if (gc->nup[i-1][0] <= 0) gc->nup[i-1][0] = 1e-20;
    if (gc->ndown[i][0] <= 0) gc->ndown[i][0] = 1e-20;
    gc->nup[i-1][1] *= gc->nup[i-1][0];
    gc->ndown[i][1] *= gc->ndown[i][0];

    if (gc->type[i] != GCX_PURE) continue;

    /* get pure state data */
    p += next;
    if (ver >= 4) { /* ring data */
      if ( 2 != sscanf(p, "%lf%lf%n",
                &gc->ring[n][0], &gc->ring[n][1], &next) ) {
        fprintf(stderr, "%s(%d): no ring data (%d) on line %d\n%s", fn, tag, n, i, s);
        break;
      }
      p += next;
      if (gc->ring[n][0] <= 0) gc->ring[n][0] = 1e-20;
      gc->ring[n][1] *= gc->ring[n][0];
    }

    if (7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf",
          &gc->nedg[n][0], &gc->ncsp[n][0], &gc->fbsm[n][0],
          &gc->nedg[n][1], &gc->ncsp[n][1],
          &gc->fbsm[n][2], &gc->fbsm[n][1]) ) {
      fprintf(stderr, "%s(%d) corrupted on line %d\n%s", fn, tag, i, s);
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
  /* should be gc->nens if successful */
  return i;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int gc_loadZrr(gc_t *gc, const char *fn, int loaddata)
{
  FILE *fp;
  int i = -1, nens, d = 0, n, ens0 = 0, n0 = 1, m = 0, ver = 0, next = 0;
  char s[512], fndef[64], sver[16];
  double tot;

  mkfnZrrdef(fn, fndef, D, gc->m - 1, gc->nmax, inode);
  xfopen(fp, fn, "r", return -1);

  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if ( s[0] != '#' || s[1] != (nmvtype ? 'R' : 'S')
    || 7 != sscanf(s + 2, "%d%d%d%d%d%d%8s%n",
                   &d, &nens, &ens0, &n0, &n, &m, sver, &next)
    || d != D || m != gc->m || ens0 != gc->ens0 ) {
    fprintf(stderr, "%s: D %d vs. %d, m %d vs %d, ens0 %d vs %d\n%s",
        fn, d, D, m, gc->m, ens0, gc->ens0, s);
    fclose(fp);
    return -1;
  }
  ver = atoi(sver[0] == 'V' ? sver + 1 : sver); /* strip V */
  if ( loaddata && (ver <= 0 || strstr(s, "dirty") || strstr(s, "DIRTY")) ) {
    fprintf(stderr, "%s: %s has no restartable data\n", fn, sver);
    loaddata = 0;
  }
  i = gc_fscanZrr(gc, fp, loaddata, fn, -1, ver, &tot);
  fclose(fp);
  printf("%4d: loaded Zrr %s, V%d, tot %g, data %d\n",
      inode, fn, ver, tot, loaddata);
  return i;
}



/* single-vertex move with a distance restraint */
INLINE int bcrstep(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xi, real r2ij[][DG_NMAX], real r2i[],
    real rc, real amp, int gauss)
{
  int i, j, n = g->n, gdirty;
  dgvs_t vs;

  DISPRNDI(i, n, x, xi, amp, gauss);
  /* check if the distance constraint is satisfied */
  if (i == i0) { /* check xi and x[j0] */
    if ( rvn_dist2(xi, x[j0]) >= rc * rc )
      return 0;
  } else if (i == j0) { /* check xi and x[i0] */
    if ( rvn_dist2(xi, x[i0]) >= rc * rc )
      return 0;
  }
  UPDGRAPHR2(i, n, g, ng, x, xi, 1, r2i, gdirty);
  if ( !gdirty ) { /* the graph is unchanged */
    rvn_copy(x[i], xi);
    UPDR2(r2ij, r2i, n, i, j);
    return 1;
  /* the new graph needs to be biconnected and without
   * the articulation pair (i0, j0) */
  } else if ( dg_biconnected(ng)
    && dg_connectedvs(ng, dgvs_mkinvset2(vs, n, i0, j0)) ) {

    dg_copy(g, ng);
    rvn_copy(x[i], xi);
    UPDR2(r2ij, r2i, n, i, j);
    return 1;
  }
  return 0;
}



#define changenapair(i, j, g, r2ij, rc) \
    changenapair_metro(i, j, g, r2ij, rc)

#if 0
/* change the distance-restrained pair */
INLINE int changenapair_simple(int *i, int *j, const dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int tmp;

  (void) g; (void) r2ij; (void) rc;
  tmp = *i, *i = *j, *j = tmp;
  return 1;
}
#endif


/* change the distance-restrained pair */
INLINE int changenapair_metro(int *i, int *j, const dg_t *g,
    real r2ij[][DG_NMAX], real rc)
{
  int ni, nj;
  dgvs_t vs;

  ni = randpair(g->n, &nj);
  if ( !(ni == *i && nj == *j) && !(ni == *j && nj == *i)
    && r2ij[ni][nj] < rc * rc
    && dg_connectedvs(g, dgvs_mkinvset2(vs, g->n, ni, nj)) ) {
    *i = ni;
    *j = nj;
    return 1;
  }
  ni = *i, *i = *j, *j = ni;
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

  *i0 = (n == 1) ? 0 : (int) (rnd0() * n); /* randomly choose a vertex as the root */
  /* xn = x0 + rc * rndball */
  rvn_inc( rvn_rndball(x[n], rc), x[*i0] );
  /* biconnectivity means to be linked to two vertices */
  for (i = 0; i < n; i++)
    if ( (r2ij[n][i] = r2ij[i][n] = rvn_dist2(x[i], x[n])) < 1 )
      deg++;
  if ( n == 1 ) deg *= 2;
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
  dgvs_t vs;

  if (Zr < 1 && rnd0() >= Zr) return 0;
  i = (rnd0() > 0.5) ? i0 : j0;
  if ( n > 2 && !dg_biconnectedvs(g, dgvs_mkinvset(vs, n, i)) )
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



/* scale the distance between x[i0] and x[j0] and test the biconnectivity */
INLINE int nmove_scale(int i0, int j0, dg_t *g, dg_t *ng,
    rvn_t *x, real *xj0, real r2ij[][DG_NMAX], real *r2i, real sr,
    int extra, double aved) /* correction for pair generation */
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
  /* if moving up to a pure state, we correct the pair generation
   * probability caused by the heat-bath method
   * the `extra' is set only in this case */
  if ( extra && rnd0() >= 1.*aved/dg_nedges(ng) ) return 0;
  dg_copy(g, ng);
  rvn_copy(x[j0], xj0);
  UPDR2(r2ij, r2i, n, j0, k);
  return 1;
}



/* check if everything is okay */
static int check(const dg_t *g, rvn_t *x, real (*r2ij)[DG_NMAX],
    const gc_t *gc, int iens0, int iens, int i0, int j0,
    double t, const char *msg)
{
  int i, j, n = g->n, err = 0, type = gc->type[iens];
  real r2, rc2 = gc->rc[iens] * gc->rc[iens];
  dgvs_t vs;

  /* check the pair distances */
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      r2 = rvn_dist2(x[i], x[j]);
      if (fabs(r2ij[i][j] - r2ij[j][i]) > 1e-6)  {
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
    if ( !dg_connectedvs(g, dgvs_mkinvset2(vs, n, i0, j0)) ) {
      fprintf(stderr, "the pair (%d, %d) is articulated\n", i0, j0);
      err++;
    }
  }
  die_if (err, "%s: t %g, iens %d (%d) -> %d (%d), n %d, err %d\n",
      msg, t, iens0, gc->type[iens0], iens, gc->type[iens], gc->n[i], err);
  return err;
}



static int mcgcr(int nmin, int nmax, int mtiers, double nsteps,
    real mcamp, int neql, double nequil, int nstsave)
{
  double t, ctot = 0, cacc = 0, prtot = 0, pracc = 0;
  int ninit, iens0, iens, ensmax, ensmin, acc, it, ieql, ned = 1, loaddata;
  int pi = 0, pj = 1; /* indices of the distance-restrained pair */
  gc_t *gc;
  dg_t *g, *ng, *g1, *ng1;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  real r2ij[DG_NMAX][DG_NMAX] = {{0}}, r2i[DG_NMAX] = {0};
  dgvs_t vs;
  const char *smove = "undefined";

  gc = gc_open(nmin, nmax, mtiers, rc0, sr0);
  loaddata = restart && !bsim0;
  gc_loadZrr(gc, fninp, loaddata);
#ifdef MPI
  /* make sure all nodes are using the same parameters */
  MPI_Bcast(gc->Zr, gc->nens, MPI_DOUBLE, MASTER, comm);
  MPI_Bcast(gc->rc, gc->nens, MPI_MYREAL, MASTER, comm);
  MPI_Bcast(gc->sr, gc->nens, MPI_MYREAL, MASTER, comm);
  MPI_Barrier(comm);
  if (inode == MASTER)
    gc_print(gc, 1, NULL, NULL, NULL);
  MPI_Barrier(comm);
#else /* OpenMP */
  #pragma omp barrier
  /* in case of OpenMP, we synchronize Zr, rc, sr manually */
  if (inode == MASTER) {
    gc_print(gc, 1, NULL, NULL, NULL);
  } else { /* nonmaster node: copy from the master */
    int i;
    for (i = 0; i < gc->nens; i++) {
      gc->Zr[i] = gc_Zr[i];
      gc->rc[i] = gc_rc[i];
      gc->sr[i] = gc_sr[i];
    }
  }
#endif

  /* we limit the initial # of vertices to 8 for larger fully-connected
   * diagrams are harder to equilibrate */
  if (nmin < 8) {
    int nn = (8 > nmax) ? nmax : 8;
    ninit = nmin + (int) (rnd0() * (nn - nmin));
  } else
    ninit = nmin;
  ensmin = nmin * mtiers;
  ensmax = nmax * mtiers;
  iens = ninit * mtiers; /* lowest ensemble */

  g = dg_open(nmax);
  ng = dg_open(nmax);
  g1 = dg_open(nmax);
  ng1 = dg_open(nmax);
  initxring(x, ninit);
  calcr2ij(r2ij, x, ninit);
  mkgraphr2ij(g, r2ij, 1, ninit);
  die_if ( !dg_biconnected(g), "n %d: not biconnected", g->n);
  printf("%4d/%-4d: n %d [%d, %d], iens %d [%d, %d]\n", inode, nnodes, g->n, nmin, nmax, iens, ensmin, ensmax);

  ieql = (neql > 0) ? 1 : neql; /* stage of equilibrations if > 0 */
  /* ieql < 0 means no data collection */

  /* main loop */
  for (it = 1, t = 1; ieql || t <= nsteps; t += 1, it++) {
    die_if (g->n < nmin || g->n > nmax, "bad n %d, t %g, iens %d, mtiers %d, ensmax %d\n", g->n, t, iens, mtiers, ensmax);
    iens0 = iens;
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
          die_if (g->n < nmin || g->n > nmax, "p2r bad n %d, t %g, iens %d\n", g->n, t, iens);
        } else if (gc->type[iens + 1] == GCX_PURE) { /* restrained up to pure */
          /*  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n  (I)
           *    -->
           *  bc(R^n) dR^n (III) */
          /* In case of heat-bath move, there is an implicit intermediate state
           *  bc(R^n) step(1 - R_{i0, j0}) c(R^n\{i0, j0}) dR^n (II)
           * The (I) --> (II) step follows from the regular scaling code.
           * If the above scaling step is accepted, then the (II) --> (III) step
           *  amounts to the following test (to remove the heat-bath bias)
           *    acc = ( rnd0() < 1. / (gc->Zr[iens + 1] * dg_nedges(ng)) );
           * where `ng' is the diagram correspond to the scaled coordinates R
           * and `1/Zr[iens + 1]' should be roughly equal to the inverse of the
           * number of edges in `ng', which is roughly n.
           * Note: if the normal Metropolis move is used for vertex removal,
           * 1/Zr[iens + 1] should be 1, which iens + 1 being a pure state */
          if ( nmvtype == 0
            || r2ij[pi][pj] * dblsqr(gc->sr[iens]) < 1 ) //dblsqr(gc->rc[iens + 1])
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i, gc->sr[iens],
               nmvtype == 1, 1. / gc->Zr[iens + 1]);
          /* the last two parameters invoke the correction of pair generation */
        } else { /* restrained up to a less restrained state */
          /*  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           *    -->
           *  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           * where rc < Rc */
          if ( r2ij[pi][pj] * dblsqr(gc->sr[iens]) < dblsqr(gc->rc[iens + 1]) )
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i, gc->sr[iens],
               0, 1.);
        }
        if (ieql >= 0) {
          gc->nup[iens][0] += 1;
          gc->nup[iens][1] += acc;
        }
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
        } else if (gc->type[iens] == GCX_PURE) { /* pure down to restrained */
          /*  bc(R^n) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n */
          /* in case of the heat-bath move, there is an intermediate state
           *  bc(R^n) step(1 - R_{i0, j0}) c(R^n\{i0, j0}) dR^n */
          if (nmvtype == 1) /* select an edge */
            ned = dg_randedge(g, &pi, &pj);
          else /* select a random pair */
            pj = randpair(g->n, &pi);

          if ( ( g->n <= 3 || dg_connectedvs(g, dgvs_mkinvset2(vs, g->n, pi, pj) ) )
            && r2ij[pj][pi] < dblsqr(gc->rc[iens - 1] * gc->sr[iens - 1])
            && (nmvtype == 0 || rnd0() < gc->Zr[iens] * ned) ) {
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i, 1. / gc->sr[iens - 1],
               0, 1.);
          }
        } else { /* restrained down to a more restrained state */
          /*  bc(R^n) step(Rc - R_{i0, j0}) c(R^n\{i0, j0}) dR^n
           *    -->
           *  bc(r^n) step(rc - r_{i0, j0}) c(r^n\{i0, j0}) dr^n
           * where Rc > rc */
          if ( r2ij[pj][pi] < dblsqr(gc->rc[iens - 1] * gc->sr[iens - 1]) )
            acc = nmove_scale(pi, pj, g, ng, x, xi, r2ij, r2i, 1. / gc->sr[iens - 1],
               0, 1.);
        }
        if (ieql >= 0) {
          gc->ndown[iens][0] += 1;
          gc->ndown[iens][1] += acc;
        }
        if ( acc ) iens--;
      }
      smove = "nmove";
    }
    else /* configuration sampling */
    {
      ctot += 1;
      if (g->n == 1) {
        cacc += 1;
      } else if (iens % mtiers == 0) { /* pure state sampling */
        int i;
        BCSTEPR2(acc, i, g->n, g, ng, x, xi, r2ij, r2i, mcamp, gaussdisp);
        cacc += acc;
      } else { /* distance-restrained sampling */
        cacc += bcrstep(pi, pj, g, ng, x, xi, r2ij, r2i,
            gc->rc[iens], mcamp, gaussdisp);
        prtot += 1;
        pracc += changenapair(&pi, &pj, g, r2ij, gc->rc[iens]);
      }
      smove = "cfg";
    }
STEP_END:
    if (nstcheck > 0 && it % nstcheck == 0)
      check(g, x, r2ij, gc, iens0, iens, pi, pj, t, smove);
    /* center the structure */
    if (it % nstcom == 0) shiftr2ij(g, x, r2ij);

    if (ieql >= 0) {
      gc->hist[iens] += 1;
      if (gc->type[iens] == GCX_PURE && rnd0() * nsted < 1)
        gc_accumdata(gc, g, t, nstcs, nstfb);
    }
    if ( ieql != 0 ) { /* equilibration */
      if ((int) fmod(t + .5, nequil) == 0) {
        if ( ieql > 0 ) { /* update parameters */
          gc_update(gc, mindata, updrc, gc->Zr, gc->rc, gc->sr);
          gc_computeZ(gc, gc->Zr, gc->rc, gc->sr);
          gc_saveZrr(gc, fnZrrtmp,  gc->Zr, gc->rc, gc->sr, 1);
          gc_print(gc, 0, NULL, NULL, NULL);
          printf("D %d, equilibration stage %d/%d, t %g\n", D, ieql, neql, t);
          gc_cleardata(gc);
          if (++ieql > neql) /* stop equilibration */
            ieql = 0;
        } else {
          printf("%d, D %d, t %g: equilibrated, n %d, iens %d\n", inode, D, t, g->n, iens);
          ieql = 0;
        }
        t = 0; it = 0; /* reset time */
      }
    } else { /* production */
      if (it % nstsave == 0 || t > nsteps - .5) {
        /* compute node-specific Z */
        gc_update(gc, mindata, updrc, gc->Zr1, gc->rc1, gc->sr1);
        gc_computeZ(gc, gc->Zr1, gc->rc1, gc->sr1);
        if (restart) {
          gc_saveZrr(gc, fnZrr, gc->Zr, gc->rc, gc->sr, 0);
          if (inode == MASTER)
            gc_saveZrr(gc, fnZrrtmp, gc->Zr1, gc->rc1, gc->sr1, 1);
        } else {
          gc_saveZrr(gc, fnZrr, gc->Zr1, gc->rc1, gc->sr1, 1);
        }
        if (inode == MASTER) {
          gc_saveZr(gc, fnZr);
          mtsave(NULL);
        }
        it = 0;
      }
    }
  } /* main loop ends here */

  printf("cacc %g, pracc %g\n", 1.*cacc/ctot, 1.*pracc/prtot);
  if (inode == MASTER)
    gc_print(gc, 0, gc->Zr1, gc->rc1, gc->sr1);
  gc_close(gc);
  dg_close(g);
  dg_close(ng);
  dg_close(g1);
  dg_close(ng1);
  return 0;
}



int main(int argc, char **argv)
{
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnodes);
#endif

  /* In MPI case, every node calls this function
   * In OpenMP case, only the master calls this function */
  doargs(argc, argv);

#ifdef _OPENMP
  nnodes = omp_get_num_threads();
  if (nnodes == 1) {
    nnodes = omp_get_max_threads();
    omp_set_num_threads(nnodes);
  }
  printf("set up %d OpenMP threads\n", nnodes);

#pragma omp parallel
  {
    inode = omp_get_thread_num();
#endif
    mtscramble(inode * 2038074743u + nnodes * time(NULL) + rngseed);
    appendfns(fninp0, fnZrr0, fnZr0);

    mcgcr(nmin, nmax, mtiers, nsteps, mcamp, neql, nequil, nstsave);

    DG_FREEMEMORIES();
    mtclosedef();
#ifdef _OPENMP
  }
#endif

#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

