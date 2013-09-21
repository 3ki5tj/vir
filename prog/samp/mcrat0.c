/* compute a virial coefficient by the ratio method
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 -march=native mcrat0.c -lm */

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




int n = 7; /* order */
double nequil = 1000000; /* number of equilibration steps */
double nsteps = 10000000;
real mcamp = 1.5;
int gdisp = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstcom = 100; /* frequency of the n move */
int nstrep = 1000000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
int cachesize = 5; /* number of recently visited diagrams */
double ratcr = 0; /* rate of coordinates replacement */
double r2cr = 0; /* variance of replaced coordinates */
char *fnout = NULL;

int bsim0 = 0;

double Bring = 1; /* to be loaded from file */
char *fnBring = NULL;



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp, "MC amplitude for system 0 (biconnected diagrams)");
  argopt_add(ao, "-L", "%d", &lookup, "lookup table (-1: automatic)");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-G", "%b", &gdisp, "normally-distributed displacement");
  argopt_add(ao, "-C", "%d", &cachesize, "cache size of recently visited diagrams");
  argopt_add(ao, "-H", "%lf", &ratcr, "rate of coordinates replacement");
  argopt_add(ao, "-U", "%lf", &r2cr, "squared radius of replaced coordinates");
  argopt_add(ao, "-o", NULL, &fnout, "output file");
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
  mcamp /= D;

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
    printf("D %d, n %d, %g steps, amp %g, nstfb %d, %d-bit, "
      "%s, %s disp, Bring %g\n",
      D, n, 1.*nsteps, mcamp, nstfb,
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



#define mkfndef(fn, fndef, d, n, inode) if (fn == NULL) { \
  sprintf(fndef, "mrD%dn%d.dat%d", d, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* save result */
static int save(const char *fn, double vir,
    const av0_t *fbsm, const av0_t *nrsm)
{
  FILE *fp;
  char fndef[64];

  mkfndef(fn, fndef, D, n, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#0 %d %d V0\n", D, n);
  fprintf(fp, "%16.0f %+17.14f %16.0f %+18.14f %+20.14e\n",
      fbsm->s, av0_getave(fbsm), nrsm->s, av0_getave(nrsm), vir);
  fclose(fp);
  printf("%4d: saved to %s, sum %g\n", inode, fn, fbsm->s);
  return 0;
}



/* load previous results */
static int load(const char *fn, av0_t *fbsm, av0_t *nrsm)
{
  FILE *fp;
  char fndef[64];
  char s[512];
  int d = 0, n1 = 0;

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
  if ( 4 != sscanf(s, "%lf%lf%lf%lf",
        &fbsm->s, &fbsm->sx, &nrsm->s, &nrsm->sx) ) {
    fprintf(stderr, "%s: corrupted\n", fn);
  }
  fbsm->sx *= fbsm->s;
  nrsm->sx *= nrsm->s;
  fclose(fp);
  printf("loaded from %s: sum %g\n", fn, fbsm->s);
  return 0;
}



/* report result */
static void report(const av0_t *fbsm, const av0_t *nrsm,
    const av0_t *cacc, const av0_t *racc,
    double t, int gmapid, const int *neval)
{
  double vir, fb, nr = av0_getave(nrsm), rv = 1;

  fb = av0_getave(fbsm);
  if (nr > 0) rv = Bring / nr;
  vir = -fb * rv;
  /* print readable result */
  printf("%d: D %d, n %d, t %g, vir %+.6e fb %+.7f nr %.6f ",
      inode, D, n, t, vir, fb, nr);
  printf("cacc %.3f, racc %.6f, ",
      av0_getave(cacc), av0_getave(racc));
  if (gmapid >= 0) printf("%d ", (int) gmapid);
  if (neval) printf("%d ", *neval);
  printf("\n");
  save(fnout, vir, fbsm, nrsm);
  if (inode == MASTER) mtsave(NULL);
}



/* compute a virial coefficient with a lookup table */
static void mcrat_lookup(int n, double nequil, double nsteps,
    real amp, int gdisp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], nx[DG_NMAX], xi;
  int i, j, pid, it, acc, nbc;
  dg_t *g, *ng;
  code_t code, ncode;
  double t, fb, nr;
  av0_t fbsm, nrsm, cacc, racc;
  unqid_t gmapid;

  die_if (n > DGMAP_NMAX, "no diagram map for n %d\n", n);

  av0_clear(&nrsm);
  av0_clear(&fbsm);
  av0_clear(&cacc);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, &fbsm, &nrsm);
  if (nnodes > 1) /* scramble the random number generator */
    mtscramble(inode * 2034091783u + time(NULL));

  g = dg_open(n);
  ng = dg_open(n);
  initx(x, n);
  mkgraph(g, x, n);
  die_if (!dg_biconnected(g),
      "%d: initial diagram not biconnected D %d\n", inode, D);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, n, g, ng, x, xi, amp, gdisp);
  printf("%d: equilibrated at t %g, nedges %d, nstfb %d\n",
      inode, nequil, dg_nedges(g), nstfb);
  dg_encode(g, &code); /* initialize the code */
  /* call _lookup function first, to initialize the table */
  nbc = dg_biconnected_lookup(g);
  gmapid = dgmap_[n].map[code];
  //nbc = dg_biconnected_lookuplow(n, gmapid);
  fb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
  nr = nbc ? dg_nring_lookuplow(n, gmapid) : 0;

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    if ( rnd0() < ratcr ) { /* particle replacement */
      /* it usually makes the simulation slower the result less precise
       * we keep it only for a psychological backup */
      dg_decode(g, &code);
      /* regenerate a configuration */
      acc = grepl0(x, nx, g, ng, r2cr);
      if ( acc ) {
        dg_encode(g, &code);
        fb = dg_hsfb_lookup(g);
        nr = dg_nring_lookup(g);
      }
      av0_add(&racc, acc);
    } else {
      /* randomly displace a particle */
      i = (int) (rnd0() * n);
      if (gdisp)
        rvn_granddisp(xi, x[i], amp);
      else
        rvn_rnddisp(xi, x[i], amp);

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
        acc = dg_biconnected_lookuplow(n, gmapid);
        if ( acc ) { /* accept the move */
          rvn_copy(x[i], xi);
          code = ncode;
          fb = dg_hsfb_lookuplow(n, gmapid);
          nr = dg_nring_lookuplow(n, gmapid);
        }
      }

      /* the cost of accumulating acc is negligible */
      av0_add(&cacc, acc);

#ifdef CHECK
      /* verify if fb has been correctly computed */
      {
        double fb1 = dg_hsfb(g), nr1 = dg_nring(g);
        dg_decode(g, &code);
        if (fabs(fb - fb1) > 1e-3 || fabs(nr - nr1) > 1e-3) {
          printf("%d: t %g, fb %d vs %d, nr %g vs %g\n",
              inode, t, fb, fb1, nr, nr1);
          dg_print(g);
          exit(1);
        }
      }
#endif
    }

    av0_add(&nrsm, nr);
    av0_add(&fbsm, fb);

    if (it % nstcom) rvn_rmcom(x, n);

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(&fbsm, &nrsm, &cacc, &racc, t, gmapid, NULL);
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



/* compute a virial coefficient by the ratio method
 * direct version without lookup table */
INLINE void mcrat_direct(int n, double nequil, double nsteps,
    real amp, int gdisp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], nx[DG_NMAX], xi;
  int i, it, nz, nnz, nbc, acc;
  double fb, nfb = 0;
  double t, nr, nnr;
  int hasnnr, hasnr;
  dg_t *g, *ng;
  av0_t fbsm, nrsm, cacc, racc;
  int ifb; /* if fb is computed in this step */
  int hasfb; /* if fb has been computed for the step */
  int hasnfb = 0; /* if fb has been computed for the MC trial */
  int neval = 0;
  dgque_t *que;

  av0_clear(&nrsm);
  av0_clear(&fbsm);
  av0_clear(&cacc);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, &fbsm, &nrsm);
  if (nnodes > 1) /* scramble the random number generator */
    mtscramble(inode * 2034091783u + time(NULL));

  ng = dg_open(n);
  g = dg_open(n);
  initx(x, n);
  mkgraph(g, x, n);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, n, g, ng, x, xi, amp, gdisp);
  printf("%d: equilibrated at t %g, nedges %d, nstfb %d\n",
      inode, nequil, dg_nedges(g), nstfb);
  die_if (!dg_biconnected(g),
      "%d, initial diagram not biconnected D %d\n", inode, D);
  nz = (dg_cliquesep(g) == 0);
  que = dgque_open(cachesize, (n * (n - 1) / 2 + DG_NMAX - 1) / DG_NMAX);
  fb = nz ? dg_hsfb_que(que, g, 1, NULL, NULL, NULL) : 0;
  hasfb = 1;
  ifb = 0;
  nr = dg_nring(g);
  hasnr = 1;

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
#ifdef CHECK
    {
      double fb1 = dg_hsfb(g), nr1 = dg_nring(g);
      int err = (fabs(nz - nz1) > 1e-3);
      if ( hasfb && fabs(fb - fb1) > 1e-3 ) err = 1;
      if ( hasnr && fabs(nr - nr1) > 1e-3 ) err = 1;
      if (err) {
        printf("%d: PRE1 t %g, hasfb %d, fb %d (%d)\n",
            inode, t, hasfb, fb, fb1);
        dg_print(g);
        exit(1);
      }
    }
#endif

    if (rnd0() < ratcr) { /* replace all coordinates */
      /* it usually makes the simulation slower the result less precise
       * we keep it only for a psychological backup */
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
      DISPRNDI(i, n, x, xi, amp, gdisp);
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
          hasnfb = 0;
          nfb = 0;
        } /* end of the general case */
        hasnnr = 0;
        nnr = 0;
      } /* end if biconnected */

      acc = nbc;
      if ( acc ) { /* accept the move */
        rvn_copy(x[i], xi);
        dg_copy(g, ng);
        nz = nnz;
        hasfb = hasnfb;
        if (hasfb) fb = nfb;
        hasnr = hasnnr;
        if (hasnr) nr = nnr;
      } /* if rejected, maintain the old `hasfb' value */

      av0_add(&cacc, acc);
    }

    /* compute fb only when necessary
     * this usually happens after an accepted move */
    ifb = (nstfb == 1 || it % nstfb == 0);
    if ( !hasfb && ifb ) {
      fb = nz ? dg_hsfb_que(que, g, 1, NULL, NULL, &neval) : 0;
      hasfb = 1;
    }

    /* compute nr only when necessary */
    if ( !hasnr && ifb ) {
      nr = dg_nring(g);
      hasnr = 1;
    }

#ifdef CHECK
    /* check if fb has been correctly computed */
    if (hasfb) {
    //if ((ifb || sys >= 2)) {
      int nz1 = (dg_cliquesep(g) == 0), ncl = dg_ncsep(g);
      int fb1 = dg_hsfb(g);
      if (!hasfb || fabs(nz - nz1) > 1e-3 || fabs(fb - fb1) > 1e-3) {
        printf("%d: t %g, acc %d, ifb %d, hasfb %d, nz %d (%d), fb %d (%d)\n",
            inode, t, acc, ifb, hasfb, nz, nz1, fb, fb1);
        dg_print(g);
        exit(1);
      }
    }
    if (hasnr) {
      double nr1 = dg_nring(g);
      if (fabs(nr - nr1) > 1e-3) {
        printf("%d: t %g, acc %d, hasnr %d, nr %g (%g)\n",
            inode, t, acc, hasnr, nr, nr1);
        dg_print(g);
        exit(1);
      }
    }
#endif

    /* accumulate data */
    if ( ifb ) {
      av0_add(&fbsm, fb);
      av0_add(&nrsm, nr);
    }

    if (it % nstcom == 0)
      rvn_rmcom(x, n); /* remove the origin to the center of mass */

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(&fbsm, &nrsm, &cacc, &racc, t, -1, &neval);
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

