/* compute a virial coefficient by the ratio method
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 -march=native mcrat0.c -lm */

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#define ZCOM_AV
#include "zcom.h"
#include "dgmapl.h" /* larger lookup table */
#include "dghash.h"
/* dg.h, ..., dgring.h are included in dgmapl.h or dghash.h if they exist
 * the preferred DG_WORDBITS are set there, if neither is available, we include
 * dgring.h directly, which in turn includes dgsc.h, dgrjw.h, etc. */
#include "dgring.h"
#include "mcutil.h"


#ifdef N
int n = N;
#else
int n = 7; /* order */
#endif
double nequil = 1000000; /* number of equilibration steps */
double nsteps = 10000000;
real mcamp = 1.5;
int gdisp = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstcom = 10000; /* frequency of the n move */
int nstrep = 1000000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
double ratcr = 0; /* rate of coordinates replacement */
double r2cr = 0; /* variance of replaced coordinates */
char *fnout = NULL;
char *prefix0 = NULL, prefix[FILENAME_MAX] = "";

int bsim0 = 0;

double Bring = 1; /* to be loaded from file */
char *fnBring = NULL;

double Zn = 0; /* partition function */

int nthreads = 1;
unsigned rngseed = 0;

int mapl_on = -1;
int mapl_kdepth = 0; /* number of links to search in the larger lookup table */

int hash_on = -1; /* 1: on, 0: off, -1: default */
int hash_bits = 0;
unsigned hash_blksz = 0;
int hash_blkmem = 0;
double hash_memmax = 0;
int hash_initls = 0; /* TODO: currently turn off initls by default mainly for compatibility
                      * Turning on initls improves the performance of memory allocation
                      * but this assumes the hash_bits and other parameters are set properly
                      * So it is safer to turn it off by default */
int auto_level = -1;
int hash_isoenum = -1;
int hash_isomax = 0; /* default value */
int hash_nocsep = 0; /* test clique separators before passing it to the hash table */

int dostat = 0; /* applies to mapl and hash */



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
  argopt_add(ao, "-H", "%lf", &ratcr, "rate of coordinates replacement");
  argopt_add(ao, "-U", "%lf", &r2cr, "squared radius of replaced coordinates");
  argopt_add(ao, "-o", NULL, &fnout, "output file");
  argopt_add(ao, "-P", NULL, &prefix0, "directory to save data");
  argopt_add(ao, "-B", "%b", &bsim0, "discard data in previous simulations");
  argopt_add(ao, "-V", "%lf", &Bring, "value of the ring integral");
  argopt_add(ao, "-I", NULL, &fnBring, "name of the virial series file");
  argopt_add(ao, "-Z", "%lf", &Zn, "partition function");
  argopt_add(ao, "--nthr", "%d", &nthreads, "number of threads");
  argopt_add(ao, "--rng", "%u", &rngseed, "set the RNG seed offset");

  /* we allow options for mapl and hash, even if they are unavailable */
  argopt_add(ao, "--mapl-mode",   "%d", &mapl_on,      "turn on/off the larger lookup table");
  argopt_add(ao, "-k",            "%d", &mapl_kdepth,  "number of links to search in the larger lookup table");
  argopt_add(ao, "--mapl-kdepth", "%d", &mapl_kdepth,  "number of links to search in the larger lookup table");

  argopt_add(ao, "--hash-mode",   "%d",   &hash_on,     "turn on/off the hash table (default: -1), see the comment of dghash_open() for detailed explanation of the hash table");
  argopt_add(ao, "--hash-bits",   "%d",   &hash_bits,   "number of bits in the key of the hash table");
  argopt_add(ao, "--hash-blksz",  "%u",   &hash_blksz,  "number of items in a link of a hash list (bucket)");
  argopt_add(ao, "--hash-blkmem", "%d",   &hash_blkmem, "number of list links to cache each time (pass to the pooled memory allocator blkmem_new() to avoid memory fragmentation), 0: default");
  argopt_add(ao, "--hash-memmax", "%lf",  &hash_memmax, "maximal memory for the hash table");
  argopt_add(ao, "--hash-initls", "%d",   &hash_initls, "initially allocate lists");
  argopt_add(ao, "--auto-level",  "%d",   &auto_level,  "automorphism level, -1: canonical label, 0: no transformation, 1: degree sequence, 2 or 3: first automorphism in the searching tree");
  argopt_add(ao, "--hash-isoenum","%d",   &hash_isoenum,"enumerate isomorphic graphs after a new graph is found, 1: yes, 0: no, -1: default");
  argopt_add(ao, "--hash-isomax", "%d",   &hash_isomax, "maximal number of items to used in the above enumeration");
  argopt_add(ao, "--hash-nocsep", "%b",   &hash_nocsep, "only pass graphs with no clique separator to the hash table");

  argopt_add(ao, "--stat",        "%b",   &dostat,      "compute statistics");

  argopt_parse(ao, argc, argv);

#ifdef N
  die_if (n != N, "n %d mismatches the predefined N %d\n", n, N);
#endif

  /* load the partition function */
  if (Zn <= 0) Zn = loadZ(D, n, NULL);
#ifdef DG_NORING
  die_if (Zn <= 0, "must supply a valid value for the partition function %g\n", Zn);
#endif /* !defined(DG_NORING) */

  /* decide whether to use lookup table or not */
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);

  /* decide if we should use the larger lookup table */
#ifdef DGMAPL_EXISTS
  /* we prefer the larger lookup table mapl to the hash table
   * due to its simplicity
   * 1. if neither mapl or hash is set, and mapl is applicable, use mapl
   * 2. if both mapl and hash are set, use mapl */
  if (mapl_on < 0) { /* default */
    if (hash_on > 0) mapl_on = 0;
    else mapl_on = (n <= DGMAPL_NMAX);
  } else if (mapl_on > 0) { /* explicitly turned */
    mapl_on = (n <= DGMAPL_NMAX);
    hash_on = 0;
  }
  /* the hash table (at leat at n <= 9) is generally not as efficient
   * as the larger lookup table, because we have to compute
   * the canonical label */
  if (mapl_on) hash_on = 0;
#else
  mapl_on = 0;
#endif

  /* decide if we should use the hash table */
#ifdef DGHASH_EXISTS
  if (hash_on != 0) {
    hash_on = 1; /* turn on the hash table generally */
    if (D > 15 && n <= 64) /* TODO: improve this */
      hash_on = 0;
  }
#else
  hash_on = 0;
#endif

#ifndef DG_NORING
  if (hash_nocsep) {
    fprintf(stderr, "hash-nocsep has little benefit without defining DG_NORING\n");
  }
#endif

  /* interval of computing fb */
  if (nstfb <= 0) {
    if (lookup || mapl_on || hash_on) {
      nstfb = 1;
    } else { /* TODO: improve this */
      nstfb = 10;
    }
  }

  /* Monte Carlo move amplitude */
  mcamp /= D;

  //if (ratcr < 0) ratcr = 0.1/n;

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
      "%s, %s disp, Bring %g\n",
      D, n, 1.*nsteps, mcamp, nstfb,
      (int) sizeof(dgword_t) * 8, lookup ? "lookup" : "direct",
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
  sprintf(fndef, "%smrD%dn%d.dat%d", prefix, d, n, inode); \
  if (inode == MASTER) fndef[strlen(fndef) - 1] = '\0'; \
  fn = fndef; }



/* save result */
static int save(const char *fn, double vir,
    const av0_t *fbsm, const av0_t *nrsm)
{
  FILE *fp;
  char fndef[64] = "";

  mkfndef(fn, fndef, D, n, inode);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#0 %d %d V1 %20.14e %20.14e\n", D, n, Bring, Zn);
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
    av0_clear(fbsm);
    av0_clear(nrsm);
    fprintf(stderr, "%s: corrupted\n", fn);
    fclose(fp);
    return -1;
  }
  fbsm->sx *= fbsm->s;
#ifndef DG_NORING
  nrsm->sx *= nrsm->s;
#else /* clear the ring data for safe */
  av0_clear(nrsm);
#endif /* !defined(DG_NORING) */

  fclose(fp);
  printf("loaded from %s: sum %g\n", fn, fbsm->s);
  return 0;
}



/* report result */
static void report(const av0_t *fbsm, const av0_t *nrsm,
    const av0_t *cacc, const av0_t *racc,
    double t, int gmapid, const int *neval)
{
  double vir, fb, nr, rv = 1;

  fb = av0_getave(fbsm);
#ifndef DG_NORING /* estimate the virial coefficient from the ring integral */
  nr = av0_getave(nrsm);
  if (nr > 0) rv = Bring / nr;
  vir = -fb * rv;
#else /* estimate the virial coefficient from the partition function */
  nr = 0;
  /* Bn = Zn <fb> * (1 - n)/ n !
   * Bn/B2^(n-1) = 2^(n-1) Zn/Z2^(n-1) */
  {
    int k;
    for (rv = 1., k = 2; k <= n; k++)
      rv *= 2./k;
  }
  rv *= Zn * (n - 1.);
#endif /* !defined(DG_NORING) */
  vir = -fb * rv;
  /* print readable result */
  printf("%d: D %d, n %d, t %g, vir %+.6e fb %+.7f nr %.6f ",
      inode, D, n, t, vir, fb, nr);
  printf("cacc %.3f, racc %.6f, rv %g, ",
      av0_getave(cacc), av0_getave(racc), rv);
  if (gmapid >= 0) printf("%d ", (int) gmapid);
  if (neval) printf("%d ", *neval);
  printf("\n");
  save(fnout, vir, fbsm, nrsm);
  if (inode == MASTER) mtsave(NULL);
}



#ifdef CHECK
  /* define the interval of checking */
  #ifndef NSTCHECK
  #define NSTCHECK 1
  #endif /* !defined(NSTCHECK) */
#endif /* defined(CHECK) */



#ifdef DGMAP_EXISTS
/* compute a virial coefficient with a lookup table */
static void mcrat_lookup(int n, double nequil, double nsteps,
    real amp, int gdisp, int nstcom)
{
  rvn_t x[DG_NMAX], nx[DG_NMAX], xi;
  int i, j, pid, it, acc;
  dg_t *g, *ng;
  dgword_t code, ncode;
  double t, fb;
#ifndef DG_NORING
  double nr = 0;
#endif
  av0_t fbsm, nrsm, cacc, racc;
  unqid_t gmapid;

  die_if (DG_N_ > DGMAP_NMAX, "no diagram map for n %d\n", DG_N_);

  av0_clear(&nrsm);
  av0_clear(&fbsm);
  av0_clear(&cacc);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, &fbsm, &nrsm);

  g = dg_open(n);
  ng = dg_open(n);
  initxring(x, DG_N_);
  mkgraph(g, x, DG_N_);
  die_if (!dg_biconnected(g),
      "%4d: initial diagram not biconnected D %d\n", inode, D);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, DG_N_, g, ng, x, xi, amp, gdisp);
  printf("%4d: equilibrated at t %g, nedges %d, rng: %p, %#x\n",
      inode, nequil, dg_nedges(g), mr_->arr, (unsigned) mr_->arr[0]);
  /* call _lookup function first, to initialize the table */
  die_if (!dg_biconnected_lookup(g),
      "initial graph (n = %d) is not biconnected\n", DG_N_);
  dg_encode(g, &code); /* initialize the code */
  gmapid = dgmap_[DG_N_].map[code];
  fb = dg_hsfb_lookuplow(DG_N_, gmapid);
#ifndef DG_NORING
  nr = dg_nring_lookuplow(DG_N_, gmapid);
#endif /* !defined(DG_NORING) */

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    if (ratcr > 0 && (ratcr >= 1 || rnd0() < ratcr)) {
      /* particle replacement */
      /* it usually makes the simulation slower the result less precise
       * we keep it only for a psychological backup */
      dg_decode(g, &code);
      /* regenerate a configuration */
      acc = grepl0(x, nx, g, ng, r2cr);
      if ( acc ) {
        dg_encode(g, &code);
        fb = dg_hsfb_lookup(g);
#ifndef DG_NORING
        nr = dg_nring_lookup(g);
#endif /* !defined(DG_NORING) */
      }
      av0_add(&racc, acc);
    } else {
      /* randomly displace a particle */
      i = (int) (rnd0() * DG_N_);
      if (gdisp)
        rvn_granddisp(xi, x[i], amp);
      else
        rvn_rnddisp(xi, x[i], amp);

      /* directly change the connectivity bitstring of the graph */
      ncode = code;
      /* for j < i pairs */
      for (pid = i - 1, j = 0; j < i; j++, pid += DG_N_ - j - 1)
        if (rvn_dist2(xi, x[j]) >= 1)
          ncode &= ~MKBIT(pid);
        else
          ncode |= MKBIT(pid);
      /* for j > i pairs */
      for (j = i + 1, pid = DG_N_*i - j*i/2; j < DG_N_; j++, pid++)
        if (rvn_dist2(xi, x[j]) >= 1)
          ncode &= ~MKBIT(pid);
        else
          ncode |= MKBIT(pid);

      if (ncode == code) { /* topology unchanged */
        acc = 1;
        rvn_copy(x[i], xi);
      } else {
        gmapid = dgmap_[DG_N_].map[ncode];
        acc = dg_biconnected_lookuplow(DG_N_, gmapid);
        if ( acc ) { /* accept the move */
          rvn_copy(x[i], xi);
          code = ncode;
          fb = dg_hsfb_lookuplow(DG_N_, gmapid);
#ifndef DG_NORING
          nr = dg_nring_lookuplow(DG_N_, gmapid);
#endif /* !defined(DG_NORING) */
        }
      }

      /* the cost of accumulating acc is negligible */
      av0_add(&cacc, acc);

#ifdef CHECK
      /* verify if fb has been correctly computed */
      if (fmod(t, NSTCHECK) < 0.1) {
        double fb1;
        dg_decode(g, &code);
        fb1 = dg_hsfb(g);
        if (fabs(fb - fb1) > 1e-3) {
          printf("%d: t %g, fb %g vs %g\n", inode, t, fb, fb1);
          dg_print(g);
          exit(1);
        }
#ifndef DG_NORING
        double nr1;
        nr1 = dg_nring(g);
        if (fabs(nr - nr1) > 1e-3) {
          printf("%d: t %g, nr %g vs %g\n", inode, t, nr, nr1);
          dg_print(g);
          exit(1);
        }
#endif /* !defined(DG_NORING) */
      }
#endif /* defined(CHECK */
    }

    av0_add(&fbsm, fb);
#ifndef DG_NORING
    av0_add(&nrsm, nr);
#endif /* !defined(DG_NORING) */

    if (it % nstcom == 0) rvn_rmcom(x, DG_N_);

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(&fbsm, &nrsm, &cacc, &racc, t, gmapid, NULL);
    }
  }
  dg_close(g);
  dg_close(ng);
}
#endif /* defined(DGMAP_EXISTS) */



#ifdef DGHASH_EXISTS
/* conditionally use the hash table to compute fb and nr */
INLINE double hashgetfbnr(dghash_t *hash, const dg_t *g, double *nr, int *ned, int *degs)
{
#ifdef DG_NORING
  int nocsep = 0;

  if ( hash_nocsep ) { /* check if there is a clique separator */
    if ( dg_cliquesep(g) ) { /* there is a clique separator */
      *nr = 0;
      return 0;
    } else { /* there is no clique separator */
      nocsep = 1;
    }
  }
  return dghash_fbnr_lookup0(hash, g, nr, nocsep, ned, degs);
#else /* !defined(DG_NORING) */
  /* if both fb and nr are needed, call dghash_fbnr_lookup0() directly */
  return dghash_fbnr_lookup0(hash, g, nr, 0, ned, degs);
#endif /* defined(DG_NORING) */
}
#endif /* defined(DGHASH_EXISTS) */



/* compute a virial coefficient by the ratio method
 * direct version without lookup table */
INLINE void mcrat_direct(int n, double nequil, double nsteps,
    real amp, int gdisp, int nstfb, int nstcom)
{
  rvn_t x[DG_NMAX], xi;
  int i, it, acc, hasfb, gdirty = 0, neval = 0;
  int ned = 0, degs[DG_NMAX];
  double t, fb, nr = 0;
  dg_t *g, *ng;
  av0_t fbsm, nrsm, cacc, racc;

  /* both the large map and hash table are shared among threads */
#ifdef DGMAPL_EXISTS
  static dgmapl_t *mapl = NULL;
#endif
#ifdef DGHASH_EXISTS
  static dghash_t *hash = NULL;
#endif

#ifdef DGMAPL_EXISTS
  if ( mapl_on ) {
#pragma omp critical
    if ( mapl == NULL ) {
      mapl = dgmapl_open(n, mapl_kdepth);
      mapl->dostat = dostat;
    }
    //printf("%4d: mapl %p\n", inode, mapl);
  }
#endif

#ifdef DGHASH_EXISTS
  if ( hash_on ) {
#pragma omp critical
    if ( hash == NULL ) {
      size_t memmax;
      /* since hash_memmax is a double, we need to do the conversion properly */
      if (hash_memmax >= (double) ((size_t) (-1) - .5))
        memmax = (size_t) (-1);
      else
        memmax = (size_t) (hash_memmax + .5);
      hash = dghash_open(n, hash_bits,
          hash_blksz, hash_blkmem, memmax, hash_initls,
          auto_level, hash_isoenum, hash_isomax);
      hash->dostat = dostat;
    }
    //printf("%4d: hash %p\n", inode, hash);
  }
#endif

  av0_clear(&nrsm);
  av0_clear(&fbsm);
  av0_clear(&cacc);
  av0_clear(&racc);
  if (!bsim0) /* try to load previous data */
    load(fnout, &fbsm, &nrsm);

  ng = dg_open(n);
  g = dg_open(n);
  initxring(x, DG_N_);
  mkgraph(g, x, DG_N_);

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    BCSTEP(acc, i, DG_N_, g, ng, x, xi, amp, gdisp);
  printf("%4d: equilibrated at t %g, nedges %d, nstfb %d, rng %p, %#x\n",
      inode, nequil, dg_nedges(g), nstfb, mr_->arr, (unsigned) mr_->arr[0]);
  die_if (!dg_biconnected(g),
      "%d, initial diagram not biconnected D %d\n", inode, D);
  fb = dg_hsfb_mixed(g);
#ifndef DG_NORING
  nr = dg_nring_mixed(g);
#endif /* !defined(DG_NORING) */
  hasfb = 1;

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
#ifdef CHECK
    if (fmod(t, NSTCHECK) < 0.1 && hasfb) {
      double fb1 = dg_hsfb(g);
      if ( fabs(fb - fb1) > 1e-3 ) {
        printf("%d: PRE1 t %g, fb %g (%g)\n", inode, t, fb, fb1);
        dg_print(g);
        exit(1);
      }
#ifndef DG_NORING
      double nr1 = dg_nring(g);
      if ( fabs(nr - nr1) > 1e-3 ) {
        printf("%d: PRE1 t %g, nr %g (%g)\n", inode, t, nr, nr1);
        dg_print(g);
        exit(1);
      }
#endif /* !defined(DG_NORING) */
    }
#endif /* defined(CHECK) */

    if (ratcr > 0 && (ratcr >= 1 || rnd0() < ratcr)) {
      /* displace one particle attached to an edge */
      acc = verepl(x, xi, g, ng, &gdirty);
      racc.sx += acc;
      racc.s += 1;
      if ( acc && gdirty ) hasfb = 0;
    } else {  /* randomly displace a particle */
      DISPRNDI(i, DG_N_, x, xi, amp, gdisp);
      UPDGRAPH(i, DG_N_, g, ng, x, xi, gdirty);
      if ( !gdirty ) {
        acc = 1;
        rvn_copy(x[i], xi);
      } else if ( (acc = dg_biconnected(ng)) != 0 ) {
        rvn_copy(x[i], xi);
        dg_copy(g, ng);
        hasfb = 0;
      }
      cacc.sx += acc;
      cacc.s += 1;
    }

    if (it % nstfb == 0) {
      if ( !hasfb ) {
        ned = -1;
        /* detect special cases, including the ring diagram
         * but it does not detect clique separators */
        if (dg_fbnr_spec0(g, &fb, &nr, &ned, degs) != 0) {
#ifdef DGMAPL_EXISTS
          if ( mapl_on && mapl != NULL ) { /* use the larger lookup table */
            fb = dgmapl_fbnr_lookup0(mapl, g, &nr, 0, &ned, degs);
          } else
#endif /* defined(DGMAPL_EXISTS) */
          {
#ifdef DGHASH_EXISTS
            if ( hash_on && hash != NULL) {
              fb = hashgetfbnr(hash, g, &nr, &ned, degs);
            } else
#endif /* defined(DGHASH_EXISTS) */
            {
              fb = dg_hsfb_mixed0(g, 0, &ned, degs);
#ifndef DG_NORING
              nr = dg_nring_mixed0(g, &ned, degs);
#endif /* !defined(DG_NORING) */
              neval++;
            }
          }
        }
        hasfb = 1;
      }
      av0_add(&fbsm, fb);
#ifndef DG_NORING
      av0_add(&nrsm, nr);
#endif
    }

#ifdef CHECK
    if (fmod(t, NSTCHECK) < 0.1) {
      /* check if fb and nr have been correctly computed */
      if (hasfb) {
        double fb1 = dg_hsfb(g);
        if (fabs(fb - fb1) > 1e-3) {
          printf("%d: t %g, acc %d, fb %g (%g)\n", inode, t, acc, fb, fb1);
          dg_print(g);
          exit(1);
        }
#ifndef DG_NORING
        double nr1 = dg_nring(g);
        if (fabs(nr - nr1) > 1e-3) {
          printf("%d: t %g, acc %d, nr %g (%g)\n", inode, t, acc, nr, nr1);
          dg_print(g);
          exit(1);
        }
#endif /* !defined(DG_NORING) */
      }
    }
#endif /* defined(CHECK) */

    if (it % nstcom == 0)
      rvn_rmcom(x, DG_N_); /* remove the origin to the center of mass */

    if (it % nstrep == 0 || t > nsteps - .5) {
      it = 0;
      report(&fbsm, &nrsm, &cacc, &racc, t, -1, &neval);
#ifdef DGMAPL_EXISTS
      if (mapl != NULL && inode == MASTER)
        dgmapl_printstat(mapl, stderr);
#endif
#ifdef DGHASH_EXISTS
      if (hash != NULL && inode == MASTER)
        dghash_printstat(hash, stderr);
#endif
    }
  }
  dg_close(g);
  dg_close(ng);

#pragma omp barrier
#ifdef DGMAPL_EXISTS
  if (mapl != NULL && inode == MASTER)
    dgmapl_close(mapl);
#endif
#ifdef DGHASH_EXISTS
  if (hash != NULL && inode == MASTER)
    dghash_close(hash);
#endif
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
  fprintf(stderr, "checking code is enabled, interval %d\n", NSTCHECK);
#endif

#ifdef _OPENMP
  if (nthreads > 1) {
    nnodes = nthreads;
    omp_set_num_threads(nnodes);
  } else {
    nnodes = omp_get_num_threads();
    if (nnodes == 1) {
      nnodes = omp_get_max_threads();
      omp_set_num_threads(nnodes);
    }
  }
  printf("set up %d OpenMP threads\n", nnodes);

#pragma omp parallel
  {
    inode = omp_get_thread_num();
#endif
    /* scramble the random number generator */
    mtscramble(inode * 2038074743u + nnodes * time(NULL) + rngseed);

#ifdef DGMAP_EXISTS
    if (lookup) {
      mcrat_lookup(n, nequil, nsteps, mcamp, gdisp, nstcom);
    } else
#endif /* defined(DGMAP_EXISTS) */
    {
      mcrat_direct(n, nequil, nsteps, mcamp, gdisp, nstfb, nstcom);
    }

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

