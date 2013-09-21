/* compute the ratio of two successive virial coefficients
 * importance sampler coupled to a uniform sampler
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 -march=native mcrat.c -lm */

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#define ZCOM_AV
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "dgring.h"
#include "mcutil.h"



#ifdef MPI
#include "mpi.h"
#define MPI_MYREAL ( (sizeof(real) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT )
MPI_Comm comm = MPI_COMM_WORLD;
#endif
#define MASTER 0
int inode = MASTER, nnodes = 1;




#define NSYS 3

int n = 7; /* order */
double nequil = 1000000; /* number of equilibration steps */
double nsteps = 10000000;
real mcamp[NSYS] = {1.5f, 0.9f, 0.9f};
int gdisp = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstcom = 100; /* frequency of the n move */
int nstrep = 1000000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
double Z[NSYS] = {1, -1, -1}; /* inverse weights (etimated partition function) */
int cachesize = 5; /* number of recently visited diagrams */
double ratcr = 0; /* rate of coordinates replacement */
char *fnout = NULL;

int bsim0 = 0;

double Bring = 1; /* to be loaded from file */
char *fnBring = NULL;



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int i;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp[0], "MC amplitude for system 0 (biconnected diagrams)");
  argopt_add(ao, "-A", "%r", &mcamp[1], "MC amplitude for system 1 (nonzero biconnected diagrams)");
  argopt_add(ao, "--A2", "%r", &mcamp[2], "MC amplitude for system 2 (importance sampler)");
  argopt_add(ao, "-Z", "%lf", &Z[1], "relative partition function of system 1");
  argopt_add(ao, "--Z2", "%lf", &Z[2], "relative partition function of system 2");
  argopt_add(ao, "-L", "%d", &lookup, "lookup table (-1: automatic)");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-G", "%b", &gdisp, "normally-distributed displacement");
  argopt_add(ao, "-C", "%d", &cachesize, "cache size of recently visited diagrams");
  argopt_add(ao, "-R", "%lf", &ratcr, "rate of coordinates replacement");
  argopt_add(ao, "-B", "%b", &bsim0, "discard data in previous simulations");
  argopt_add(ao, "-V", "%lf", &Bring, "value of the ring integral");
  argopt_add(ao, "-I", NULL, &fnBring, "name of the virial series file");
  argopt_parse(ao, argc, argv);

  /* decide whether to use lookup table or not */
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);

  /* interval of computing fb */
  if (nstfb <= 0) {
    if (lookup) {
      nstfb = 1;
    } else { /* TODO: improve this */
      if (D >= 15) nstfb = 10;
      else nstfb = 100;
    }
  }

  /* Monte Carlo move amplitude */
  for (i = 0; i < NSYS; i++) mcamp[i] /= D;

  /* default weights */
  if (Z[1] < 0) {
    static double Z1[] = {1, 1, 1, 1, 0.29*D/3, 0.084*D/3, 0.044*D/3, 0.032*D/3, 0.026*D/3, 0.022*D/3, 0.019*D/3};
    Z[1] = (n > 10) ? 0.2*D/3/n : Z1[n];
  }
  if (Z[2] < 0) {
    static double Z2[] = {1, 1, 1, 1, 0.53, 0.3, 0.2, 0.17, 0.16, 0.158};
    Z[2] = (n > 9) ? 0.15 : Z2[n];
  }

  //if (ratcr < 0) ratcr = 0.1/n;

  if ( !argopt_isset(ao, Bring) ) {
    if (inode == 0)
      loadBring(D, n, n, &Bring, fnBring);
#ifdef MPI
    MPI_Bcast(&Bring, 1, MPI_DOUBLE, MASTER, comm);
#endif
  }

  if (inode == MASTER) {
    argopt_dump(ao);
    printf("D %d, n %d, %g steps, amp %g/%g/%g, nstfb %d, %d-bit, "
      "%s, %s disp, Bring %g\n",
      D, n, 1.*nsteps, mcamp[0], mcamp[1], mcamp[2], nstfb,
      (int) sizeof(code_t) * 8, lookup ? "lookup" : "direct",
      gdisp ? "Gaussian" : "uniform", Bring);
  } else { /* change file names for slave nodes */
    if (fnout) fnout = fnappend(fnout, inode);
  }
  argopt_close(ao);
#ifdef MPI
  MPI_Barrier(comm);
#endif
}



#define WTSYS2(wt, fb) if ((wt = fabs(fb)) < 1e-3) wt = 1;

/* change `sys' between three systems of biconnected diagrams
 * 1. all
 * 2. those without clique separators
 * 3. those with |fb| as the weight  */
INLINE int movesys(int sys, double fb, int nz, int ifb,
    const double *Z, av0_t *tacc)
{
  double r, wt;

  if (sys == 0) {
    if ( nz ) { /* only a `nz' configuration can switch */
      if (rnd0() < 0.5) { /* switch sys from 0 to 1 */
        if ((r = Z[0]/Z[1]) >= 1 || rnd0() < r)
          sys = 1;
      } else { /* switch sys from 0 to 2 */
        if ( ifb ) {
          WTSYS2(wt, fb);
          if ((r = wt*Z[0]/Z[2]) >= 1 || rnd0() < r)
            sys = 2;
        }
      }
    }
  } else if (sys == 1) {
    if (rnd0() < 0.5) { /* switch sys from 1 to 0 */
      if ((r = Z[1]/Z[0]) >= 1 || rnd0() < r)
        sys = 0;
    } else { /* switch sys from 1 to 2 */
      if ( ifb ) {
        WTSYS2(wt, fb);
        if ((r = wt*Z[1]/Z[2]) >= 1 || rnd0() < r)
          sys = 2;
        av0_add(&tacc[0], sys == 2);
      }
    }
  } else { /* sys 2 */
    if (rnd0() < 0.5) { /* switch sys from 2 to 1 */
      WTSYS2(wt, fb);
      if ((r = Z[2]/(Z[1]*wt)) >= 1 || rnd0() < r)
        sys = 1;
      av0_add(&tacc[1], sys == 1);
    } else { /* switch sys from 2 to 0 */
      WTSYS2(wt, fb);
      if ((r = Z[2]/(Z[0]*wt)) >= 1 || rnd0() < r)
        sys = 0;
    }
  }
  return sys;
}



/* change `sys' between three systems of biconnected diagrams
 * 1. all
 * 2. those without clique separators
 * 3. those with |fb| as the weight */
INLINE int movesys_heatbath(int sys, double fb,
    const double *Z, av0_t *tacc)
{
  double proba[4], r, wt;
  int i;

  /* compute the relative weights */
  proba[0] = 0;
  proba[1] = 1./Z[0]; /* system 0 accepts all diagrams */
  proba[2] = proba[1] + 1./Z[1];
  WTSYS2(wt, fb);
  proba[3] = proba[2] + wt/Z[2];

  r = rnd0() * proba[3];
  if (r > proba[2]) i = 2;
  else if (r > proba[1]) i = 1;
  else i = 0;

  /* compute the acceptance probabilities */
  if (sys == 1) {
    av0_add(&tacc[0], i == 2);
  } else if (sys == 2) {
    av0_add(&tacc[1], i == 1);
  }
  return i;
}



#define mkfndef(fn, fndef, d, n, inode) if (fn == NULL) { \
  sprintf(fndef, "mrD%dn%d.dat%d", d, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* compute the virial coefficients */
static void compute(double vir[], double fb[], double ta[],
    const av0_t *fbsm, const av0_t *nzsm,
    const av0_t *nrsm, const av0_t *tacc)
{
  double nz = av0_getave(nzsm), nr = av0_getave(nrsm), rv = 1;
  int i;

  ta[0] = av0_getave(&tacc[0]);
  ta[1] = av0_getave(&tacc[1]);
  fb[0] = av0_getave(&fbsm[0]);
  /* multiple histogram: fb = sum_k fbsm[k] / sum_k Z0/Zk hist[k] */
  fb[1] = nz * (fbsm[0].sx + fbsm[1].sx) / (nz * fbsm[0].s + fbsm[1].s);
  if (ta[1] > 0)
    fb[2] = av0_getave(&fbsm[2]) * (Z[2] / Z[1]) * (ta[0] / ta[1]) * nz;
  else
    fb[2] = 0;
  if (nr > 0) rv = Bring / nr;
  for (i = 0; i < NSYS; i++) vir[i] = -fb[i] * rv;
}



/* save result */
static int save(const char *fn, const double *vir, const av0_t *fbsm,
    const av0_t *nzsm, const av0_t *nrsm, const av0_t *tacc)
{
  FILE *fp;
  char fndef[64];

  mkfndef(fn, fndef, D, n, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#M %d %d 1 %.14e %.14e V0\n", D, n, Z[1], Z[2]);
  fprintf(fp, "%16.0f %+17.14f %16.0f %+17.14f %16.0f %17.14f "
      "%16.0f %+17.14f %16.0f %+18.14f %16.0f %16.14f %16.0f %16.14f "
      "%+20.14e %+20.14e %+20.14e\n",
      fbsm[0].s, av0_getave(&fbsm[0]),
      fbsm[1].s, av0_getave(&fbsm[1]),
      fbsm[2].s, av0_getave(&fbsm[2]),
      nzsm->s, av0_getave(nzsm), nrsm->s, av0_getave(nrsm),
      tacc[0].s, av0_getave(&tacc[0]), tacc[1].s, av0_getave(&tacc[1]),
      vir[0], vir[1], vir[2]);
  fclose(fp);
  printf("%4d: saved to %s, sum %g,%g,%g\n", inode, fn, fbsm[0].s, fbsm[1].s, fbsm[2].s);
  return 0;
}



/* load previous results */
static int load(const char *fn, av0_t *fbsm, av0_t *nzsm,
    av0_t *nrsm, av0_t *tacc)
{
  FILE *fp;
  char fndef[64];
  char s[512];
  int d, n1, next, i;

  mkfndef(fn, fndef, D, n, inode);
  xfopen(fp, fn, "r", return -1);
  if ( fgets(s, sizeof s, fp) == NULL
    || 2 != sscanf(s + 2, "%d%d", &d, &n1)
    || d != D || n1 != n ) {
    fprintf(stderr, "%s: bad info D %d vs %d, n %d vs %d\n%s",
        fn, d, D, n1, n, s);
    fclose(fp);
    return -1;
  }
  if ( fgets(s, sizeof s, fp) == NULL ) {
    fclose(fp);
    return -1;
  }
  if ( 14 != sscanf(s, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
        &fbsm[0].s, &fbsm[0].sx, &fbsm[1].s, &fbsm[1].sx,
        &fbsm[2].s, &fbsm[2].sx, &nzsm->s, &nzsm->sx, &nrsm->s, &nrsm->sx,
        &tacc[0].s, &tacc[0].sx, &tacc[1].s, &tacc[1].sx) ) {
    fprintf(stderr, "%s: corrupted\n", fn);
  }
  for (i = 0; i < NSYS; i++) fbsm[i].sx *= fbsm[i].s;
  nzsm->sx *= nzsm->s;
  nrsm->sx *= nrsm->s;
  for (i = 0; i < 2; i++) tacc[i].sx *= tacc[i].s;
  fclose(fp);
  printf("loaded from %s: sum %g, %g, %g\n", fn, fbsm[0].s, fbsm[1].s, fbsm[2].s);
  return 0;
}



/* report result */
static void report(const av0_t *fbsm, const av0_t *nzsm, const av0_t *nrsm,
    const av0_t *tacc, const av0_t *cacc, const av0_t *racc,
    double t, int gmapid, const int neval[])
{
  double vir[3], fb[3], ta[2];

  compute(vir, fb, ta, fbsm, nzsm, nrsm, tacc);
  /* print readable result */
  printf("%d: D %d, n %d, t %g, vir %+.6e,%+.6e,%+.6e fb %+.7f,%+.7f,%+.7f\n",
      inode, D, n, t, vir[0], vir[1], vir[2], fb[0], fb[1], fb[2]);
  printf("  nz %.6f, his 1,%.4f,%.4f; tacc %.6f %.6f, "
      "cacc %.3f/%.3f/%.3f, racc %.3f, ",
      av0_getave(nzsm), fbsm[1].s / fbsm[0].s,
      fbsm[2].s / fbsm[0].s, ta[0], ta[1],
      av0_getave(&cacc[0]), av0_getave(&cacc[1]),
      av0_getave(&cacc[2]), av0_getave(racc));
  if (gmapid >= 0) printf("%d ", (int) gmapid);
  if (neval) printf("%d, %d, %d ", neval[0], neval[1], neval[2]);
  printf("\n");
  save(fnout, vir, fbsm, nzsm, nrsm, tacc);
  if (inode == MASTER) mtsave(NULL);
}



/* compute the sign of the virial coefficient with a lookup table */
static void mcrat_lookup(int n, double nequil, double nsteps,
    real amp[], int gdisp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], nx[DG_NMAX], xi;
  int i, j, pid, it, nz, nnz, nbc;
  dg_t *g, *ng;
  code_t code, ncode;
  double t, nfb = 0, fb, wt, nwt, nr;
  av0_t fbsm[NSYS], nzsm, nrsm, tacc[2], cacc[NSYS], racc;
  int sys, sys0, acc;
  unqid_t gmapid;

  die_if (n > DGMAP_NMAX, "no diagram map for n %d\n", n);

  av0_clear(&nrsm);
  av0_clear(&nzsm);
  for (i = 0; i < NSYS; i++) {
    av0_clear(&fbsm[i]);
    av0_clear(&cacc[i]);
  }
  for (i = 0; i < 2; i++) av0_clear(&tacc[i]);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, fbsm, &nzsm, &nrsm, tacc);
  if (nnodes > 1) /* scramble the random number generator */
    mtscramble(inode * 2034091783u + time(NULL));

  g = dg_open(n);
  ng = dg_open(n);
  initx(x, n);
  mkgraph(g, x, n);
  die_if (!dg_biconnected(g),
      "%d: initial diagram not biconnected D %d\n", inode, D);
  fb = dg_hsfb_lookup(g);
  nr = dg_nring_lookup(g);
  nz = 1; /* if (lookup), nz means fb != 0,
             otherwise nz means no clique separator */
  sys = 0;

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, n, g, ng, x, xi, amp[0], gdisp);
  printf("%d: equilibrated at t %g, nedges %d, nstfb %d\n",
      inode, nequil, dg_nedges(g), nstfb);
  dg_encode(g, &code); /* initialize the code */
  gmapid = dgmap_[n].map[code];
  nbc = dg_biconnected_lookuplow(n, gmapid);
  fb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
  nr = nbc ? dg_nring_lookuplow(n, gmapid) : 0;
  nz = (fabs(fb) > 1e-3);

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    sys0 = sys;

    if (sys == 0 && rnd0() < ratcr) { /* particle replacement */
      dg_decode(g, &code);
      /* regenerate a configuration */
      acc = grepl(x, nx, g, ng);
      if ( acc ) {
        dg_encode(g, &code);
        fb = dg_hsfb_lookup(g);
        nz = (fabs(fb) > 1e-3);
        nr = dg_nring_lookup(g);
      }
      av0_add(&racc, acc);
    } else {
      /* randomly displace a particle */
      i = (int) (rnd0() * n);
      if (gdisp)
        rvn_granddisp(xi, x[i], amp[sys0]);
      else
        rvn_rnddisp(xi, x[i], amp[sys0]);

      /* directly change the connectivity bitstring of the graph */
      ncode = code;
      /* for j < i pairs */
      for (pid = i - 1, j = 0; j < i; j++, pid += n - j - 1)
        if (rvn_dist2(xi, x[j]) >= 1)
          ncode &= ~(1u << pid);
        else
          ncode |= (1u << pid);
      /* for j > i pairs */
      for (j = i + 1, pid = n*i - j*i/2; j < n; j++, pid++)
        if (rvn_dist2(xi, x[j]) >= 1)
          ncode &= ~(1u << pid);
        else
          ncode |= (1u << pid);

      if (ncode == code) { /* topology unchanged */
        acc = 1;
        rvn_copy(x[i], xi);
      } else {
        gmapid = dgmap_[n].map[ncode];
        nbc = dg_biconnected_lookuplow(n, gmapid);
        /* the weight `fb' is computed from the lookup table
         * for a small n, so it can be updated in every step */
        nfb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
        nnz = (fabs(nfb) > 1e-3);

        if (sys0 == 0) {
          acc = nbc;
        } else if (sys0 == 1) {
          acc = nnz;
        } else { /* sys0 == 2 */
          if (nnz) {
            WTSYS2(wt, fb);
            WTSYS2(nwt, nfb);
            if (nwt >= wt) acc = 1;
            else acc = (rnd0() < 1.*nwt/wt);
          } else acc = 0;
        }

        if ( acc ) { /* accept the move */
          rvn_copy(x[i], xi);
          code = ncode;
          nz = nnz;
          fb = nfb;
          nr = dg_nring_lookuplow(n, gmapid);
        }
      }

      /* change `sys', this costs about 16% of simulation time */
      if ( nz ) /* if fb == 0, then it only belongs to system 0 */
        sys = movesys_heatbath(sys0, fb, Z, tacc);

      /* the cost of accumulating acc is negligible */
      av0_add(&cacc[sys0], acc);

#ifdef CHECK
      /* verify if fb has been correctly computed */
      {
        int fb1, nfb1;
        double nr1, nnr1;
        dg_decode(g, &code);
        dg_decode(ng, &ncode);
        fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
        nr1 = dg_nring(g), nnr1 = dg_nring(ng);
        if (fb != fb1 || nr != nr1) {
          printf("%d: t %g, sys %d, nz %d (%d), fb %d (%d), nr %g (%g)"
              "acc %d, nfb %d (%d)\n",
              inode, t, sys, nz, fb != 0,
              fb, fb1, nr, nr1, acc, nfb, nfb1);
          dg_print(g);
          exit(1);
        }
      }
#endif
    }

    if (sys == 0) {
      av0_add(&nzsm, nz);
      av0_add(&nrsm, nr);
    }
    if (sys <= 1) {
      av0_add(&fbsm[sys], fb);
    } else {
      WTSYS2(wt, fb);
      av0_add(&fbsm[sys], fb/wt);
    }

    if (it % nstcom) rvn_rmcom(x, n);

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(fbsm, &nzsm, &nrsm, tacc, cacc, &racc, t, gmapid, NULL);
    }
  }
  dg_close(g);
  dg_close(ng);
}



typedef struct {
  struct {
    double val; /* e.g., the value of fb */
    code_t c[DG_NMAX/2]; /* code */
  } *arr;
  int nmax;
  int cnt; /* number of entries */
  int i; /* current point from 0 to cnt - 1 */
  int k; /* number of code_t */
} dgque_t;



INLINE dgque_t *dgque_open(int cnt, int k)
{
  dgque_t *q;

  xnew(q, 1);
  q->nmax = cnt;
  q->cnt = 0;
  q->k = k;
  die_if (k > DG_NMAX/2, "k %d, n %d is too large", k, n);
  xnew(q->arr, q->nmax);
  return q;
}



INLINE void dgque_close(dgque_t *q)
{
  free(q->arr);
  free(q);
}



INLINE int dgque_find(const dgque_t *q, const code_t *c, double *x)
{
  int i, j, cnt = q->cnt, k = q->k;

  for (i = 0; i < cnt; i++) {
    for (j = 0; j < k; j++)
      if (q->arr[i].c[j] != c[j])
        break;
    if (j == k) {
      *x = q->arr[i].val;
      return i;
    }
  }
  return -1;
}



INLINE void dgque_add(dgque_t *q, code_t *c, double x)
{
  int i, j, k = q->k;

  if (q->cnt < q->nmax) {
    i = q->cnt;
    q->arr[i].val = x;
    for (j = 0; j < k; j++) q->arr[i].c[j] = c[j];
    if ((q->cnt = i + 1) == q->nmax)
      q->i = 0; /* set point */
  } else { /* replace the q->i th point */
    i = q->i;
    q->arr[i].val = x;
    for (j = 0; j < k; j++) q->arr[i].c[j] = c[j];
    q->i = (i + 1) % q->nmax;
  }
}



/* compute fb with a short history list
 * nocsep: no clique separator in the graph */
INLINE double dg_hsfb_que(dgque_t *q, const dg_t *g, int nocsep,
    int *ned, int *degs, int *neval)
{
  double fb = 0;
  static code_t c[DG_NMAX/2];

  dg_encode(g, c);
  if (dgque_find(q, c, &fb) < 0) {
    /* hsfb_mixed() is much faster than hsfb_rjw() in high dimensions D
     * because when the number of edges ~ n, rhsc() is faster than hsfb() */
    fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
    dgque_add(q, c, fb);
    if (neval) (*neval)++;
  }
  return fb;
}



/* compute the sign of the virial coefficient by importance sampling
 * direct version without lookup table */
INLINE void mcrat_direct(int n, double nequil, double nsteps,
    real amp[], int gdisp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], nx[DG_NMAX], xi;
  int i, it, nz, nnz, nbc;
  double wt, nwt, fb, nfb = 0;
  double t, nr, nnr;
  int hasnnr, hasnr;
  dg_t *g, *ng;
  av0_t fbsm[NSYS], nzsm, nrsm, tacc[2], cacc[NSYS], racc;
  int sys, sys0, acc;
  int ifb; /* if fb is randomly computed */
  int hasfb; /* if fb has been computed for the step */
  int hasnfb = 0; /* if fb has been computed for the MC trial */
  int neval[3] = {0, 0, 0};
  dgque_t *que;

  av0_clear(&nrsm);
  av0_clear(&nzsm);
  for (i = 0; i < NSYS; i++) {
    av0_clear(&fbsm[i]);
    av0_clear(&cacc[i]);
  }
  for (i = 0; i < 2; i++) av0_clear(&tacc[i]);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, fbsm, &nzsm, &nrsm, tacc);
  if (nnodes > 1) /* scramble the random number generator */
    mtscramble(inode * 2034091783u + time(NULL));

  ng = dg_open(n);
  g = dg_open(n);
  initx(x, n);
  mkgraph(g, x, n);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, n, g, ng, x, xi, amp[0], gdisp);
  printf("%d: equilibrated at t %g, nedges %d, nstfb %d\n",
      inode, nequil, dg_nedges(g), nstfb);
  /* initialize state variables */
  nbc = dg_biconnected(g);
  die_if (!nbc,
      "%d, initial diagram not biconnected D %d\n", inode, D);
  nz = (dg_cliquesep(g) == 0);
  que = dgque_open(cachesize, (n * (n - 1) / 2 + DG_NMAX - 1) / DG_NMAX);
  fb = nz ? dg_hsfb_que(que, g, 1, NULL, NULL, NULL) : 0;
  hasfb = 1;
  ifb = 0;
  nr = dg_nring(g);
  hasnr = 1;
  sys = 0;

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* save the system, 0: biconnected; 1: nocsep; 2: w ~ 1/|fb| */
    sys0 = sys;

#ifdef CHECK
    {
      int nz1 = (dg_cliquesep(g) == 0), ncl = dg_ncsep(g);
      int fb1 = dg_hsfb(g);
      double nr1 = dg_nring(g);
      int err = (nz != nz1);
      if ( hasfb && fb != fb1 ) err = 1;
      if ( hasnr && nr != nr1 ) err = 1;
      if (err) {
        printf("%d: PRE1 t %g, sys %d, hasfb %d, nz %d (%d), ncl %d, fb %d (%d)\n",
            inode, t, sys, hasfb, nz, nz1, ncl, fb, fb1);
        dg_print(g);
        exit(1);
      }
    }
#endif

    if (sys0 == 0 && rnd0() < ratcr) { /* replace all coordinates */
      /* regularly regenerate a configuration */
      acc = grepl(x, nx, g, ng);
      if ( acc ) {
        nz = (dg_cliquesep(g) == 0);
        hasfb = 0;
        fb = 0;
        hasnr = 0;
        nr = 0;
      }
      av0_add(&racc, acc);
    } else {  /* randomly displace a particle */
      DISPRNDI(i, n, x, xi, amp[sys0], gdisp);
      UPDGRAPH(i, n, g, ng, x, xi);

      /* try to avoid the expensive computation of fb */
      nbc = dg_biconnected(ng);
      if ( !nbc ) {
        nnz = 0;
        nfb = 0;
        nnr = 0;
        hasnfb = 1;
        hasnnr = 1;
      } else { /* biconnected */
        int ned = -1;
        static int degs[DG_NMAX];
        int err;
        double sc;

        /* detect special cases, including ring diagrams
         * and clique separators */
        sc = dg_rhsc_spec0(ng, 0, &ned, degs, &err);
        if (err == 0) { /* special case worked */
          hasnfb = 1;
          if (fabs(sc) > 1e-3) {
            nnz = 1;
            nfb = DG_SC2FB(sc, ned);
          } else {
            nnz = 0;
            nfb = 0;
          }
        } else { /* general case */
          nnz = 1; /* no clique separator */
          /* we need to compute `nfb' if `sys' == 3 */
          if ( sys0 >= 2 ) { /* if we must compute the new fb */
            hasnfb = 1;
            /* if the connectivity is unchanged, use the old fb */
            if (hasfb && g->c[i] == ng->c[i]) {
              nfb = fb;
            } else {
              nfb = dg_hsfb_que(que, ng, 1, &ned, degs, &neval[sys0]);
            }
          } else { /* indicate fb has not been computed */
            hasnfb = 0;
            nfb = 0;
          }
        } /* end of the general case */
        hasnnr = 0;
        nnr = 0;
      } /* end if biconnected */

      if (sys0 == 0) {
        acc = nbc;
      } else if (sys0 == 1) {
        acc = nnz;
      } else { /* sys0 == 2 */
        if (nnz) {
          WTSYS2(wt, fb);
          WTSYS2(nwt, nfb);
          if (nwt >= wt) acc = 1;
          else acc = (rnd0() < 1.*nwt/wt);
        } else acc = 0;
      }

      if ( acc ) { /* accept the move */
        rvn_copy(x[i], xi);
        dg_copy(g, ng);
        nz = nnz;
        hasfb = hasnfb;
        if (hasfb) fb = nfb;
        hasnr = hasnnr;
        if (hasnr) nr = nnr;
      } /* if rejected, maintain the old `hasfb' value */

      av0_add(&cacc[sys], acc);
    }

    /* compute fb only when necessary
     * this usually happens after an accepted move */
    ifb = (nstfb == 1 || rnd0() < 1./nstfb);
    if ( !hasfb && (ifb || sys0 >= 2) ) {
      fb = nz ? dg_hsfb_que(que, g, 1, NULL, NULL, &neval[sys0]) : 0;
      hasfb = 1;
    }

    /* change `sys' */
    sys = movesys(sys0, fb, nz, ifb, Z, tacc);

    /* compute nr only when necessary */
    if ( !hasnr && ifb && sys == 0) {
      nr = dg_nring(g);
      hasnr = 1;
    }

#ifdef CHECK
    /* check if fb has been correctly computed */
    if (hasfb) {
    //if ((ifb || sys >= 2)) {
      int nz1 = (dg_cliquesep(g) == 0), ncl = dg_ncsep(g);
      int fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
      if (!hasfb || nz != nz1 || fb != fb1) {
        printf("%d: t %g, acc %d, ifb %d, sys %d/%d, hasfb %d, nz %d (%d), "
            "ncl %d, fb %d (%d), hasnfb %d, nfb %d (%d)\n",
            inode, t, acc, ifb, sys0, sys, hasfb, nz, nz1,
            ncl, fb, fb1, hasnfb, nfb, nfb1);
        dg_print(g);
        exit(1);
      }
    }
    if (hasnr) {
      double nr1 = dg_nring(g), nnr1 = dg_nring(ng);
      if (nr != nr1) {
        printf("%d: t %g, acc %d, sys %d/%d, hasnr %d, "
            "nr %g (%g), hasnnr %d, nnr %g (%g)\n",
            inode, t, acc, sys0, sys, hasnr,
            nr, nr1, hasnnr, nnr, nnr1);
        dg_print(g);
        exit(1);
      }
    }
#endif

    /* accumulate data */
    if (ifb || sys0 >= 2) {
      if (sys == 0) av0_add(&fbsm[0], nz);
      if (sys <= 1) {
        av0_add(&fbsm[sys], fb);
      } else {
        WTSYS2(wt, fb);
        av0_add(&fbsm[sys], fb/wt);
      }
    }

    if (ifb && sys == 0) av0_add(&nrsm, nr);

    if (it % nstcom == 0)
      rvn_rmcom(x, n); /* remove the origin to the center of mass */

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(fbsm, &nzsm, &nrsm, tacc, cacc, &racc, t, -1, neval);
    }
  }
  dg_close(g);
  dg_close(ng);
  dgque_close(que);
}



int main(int argc, char **argv)
{
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnodes);
#endif
  doargs(argc, argv);

#ifdef CHECK
  fprintf(stderr, "checking code is enabled\n");
#endif

  if (lookup)
    mcrat_lookup(n, nequil, nsteps, mcamp, gdisp, nstfb, nstcom);
  else
    mcrat_direct(n, nequil, nsteps, mcamp, gdisp, nstfb, nstcom);

#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

