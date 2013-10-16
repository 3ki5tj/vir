/* compute a virial coefficient by the ratio method
 * large lookup table method (due to Labik) for n = 9
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 -march=native mcrat0.c -lm */

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#define ZCOM_AV
#include "zcom.h"

#define NMAX 64 /* force to use the 64 bit representation */
#include "dgmapl.h"
#include "mcutil.h"



int n = 9; /* order */
double nequil = 1000000; /* number of equilibration steps */
double nsteps = 10000000;
real mcamp = 1.5;
int nstfb = 0; /* interval of evaluting the weight */
int nstcom = 10000; /* frequency of the n move */
int nstrep = 1000000000; /* interval of reporting */
double ratcr = 0; /* rate of coordinates replacement */
double r2cr = 0; /* variance of replaced coordinates */
char *fnout = NULL;
char *prefix0 = NULL, prefix[FILENAME_MAX] = "";
int kdepth = 0; /* number of links to search in the lookup table */

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
  argopt_add(ao, "-k", "%d", &kdepth, "number of links to search in the lookup table");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-H", "%lf", &ratcr, "rate of coordinates replacement");
  argopt_add(ao, "-U", "%lf", &r2cr, "squared radius of replaced coordinates");
  argopt_add(ao, "-o", NULL, &fnout, "output file");
  argopt_add(ao, "-P", NULL, &prefix0, "directory to save data");
  argopt_add(ao, "-B", "%b", &bsim0, "discard data in previous simulations");
  argopt_add(ao, "-V", "%lf", &Bring, "value of the ring integral");
  argopt_add(ao, "-I", NULL, &fnBring, "name of the virial series file");
  argopt_parse(ao, argc, argv);

  /* interval of computing fb */
  if (nstfb <= 0) {
    nstfb = 1;
  }

  /* Monte Carlo move amplitude */
  mcamp /= D;

  if ( !argopt_isset(ao, Bring) ) {
    if (inode == MASTER)
      loadBring(D, n, n, &Bring, fnBring);
#ifdef MPI
    MPI_Bcast(&Bring, 1, MPI_DOUBLE, MASTER, comm);
#endif
  }

  if (prefix0) sprintf(prefix, "%s/", prefix0);

  if (inode == MASTER) {
    argopt_dump(ao);
    printf("D %d, n %d, %g steps, amp %g, nstfb %d, %d-bit, "
      "Lookup, Bring %g\n",
      D, n, 1.*nsteps, mcamp, nstfb,
      (int) sizeof(code_t) * 8, Bring);
  } else { /* change file names for slave nodes */
    if (fnout) fnout = fnappend(fnout, inode);
  }
  argopt_close(ao);
#ifdef MPI
  MPI_Barrier(comm);
#endif
}



#define mkfndef(fn, fndef, d, n, inode) if (fn == NULL) { \
  sprintf(fndef, "%smrD%dn%d.dat%d", prefix, d, n, inode); \
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
  char fndef[64] = "";
  static char s[512];
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
    double t)
{
  double vir, fb, nr = av0_getave(nrsm), rv = 1;

  fb = av0_getave(fbsm);
  if (nr > 0) rv = Bring / nr;
  vir = -fb * rv;
  /* print readable result */
  printf("%d: D %d, n %d, t %g, vir %+.6e fb %+.7f nr %.6f ",
      inode, D, n, t, vir, fb, nr);
  printf("cacc %.3f, racc %.6f\n",
      av0_getave(cacc), av0_getave(racc));
  save(fnout, vir, fbsm, nrsm);
  if (inode == MASTER) mtsave(NULL);
}



/* compute a virial coefficient with a lookup table */
static void mcrat_Lookup(int n, double nequil, double nsteps,
    real amp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], xi;
  int i, it, acc, gdirty = 0, dirty = 0;
  dg_t *g, *ng;
  double t, fb, nr;
  av0_t fbsm, nrsm, cacc, racc;

  die_if (n > DGMAPL_NMAX, "no diagram map for n %d\n", n);

  av0_clear(&nrsm);
  av0_clear(&fbsm);
  av0_clear(&cacc);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, &fbsm, &nrsm);
  /* scramble the random number generator */
  mtscramble(inode * 2038074743u + nnodes * time(NULL));

  g = dg_open(n);
  ng = dg_open(n);
  initxring(x, n);
  mkgraph(g, x, n);
  die_if (!dg_biconnected(g),
      "%4d: initial diagram not biconnected D %d\n", inode, D);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, n, g, ng, x, xi, amp, 0);
  printf("%4d: equilibrated at t %g, nedges %d\n",
      inode, nequil, dg_nedges(g));
  die_if (!dg_biconnected(g),
      "initial graph (n = %d) not biconnected\n", g->n);
  fb = dg_fbnr_Lookup(g, kdepth, &nr);

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* randomly displace a particle */
    i = (int) (rnd0() * n);
    rvn_rnddisp(xi, x[i], amp);

    UPDGRAPH(i, n, g, ng, x, xi, gdirty);

    if ( !gdirty ) { /* topology unchanged */
      acc = 1;
      rvn_copy(x[i], xi);
    } else if ( (acc = dg_biconnected(ng)) ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      dirty = 1;
    }

    if (it % nstfb == 0) {
      if (dirty)
        fb = dg_fbnr_Lookup(g, kdepth, &nr);
      dirty = 0;
    }
    av0_add(&cacc, acc);
    av0_add(&nrsm, nr);
    av0_add(&fbsm, fb);

    if (it % nstcom == 0) rvn_rmcom(x, n);

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(&fbsm, &nrsm, &cacc, &racc, t);
    }
  }
  dg_close(g);
  dg_close(ng);
}



int main(int argc, char **argv)
{
#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnodes);
#endif

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

    mcrat_Lookup(n, nequil, nsteps, mcamp, nstfb, nstcom);

#ifdef _OPENMP
  }
#endif

#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

