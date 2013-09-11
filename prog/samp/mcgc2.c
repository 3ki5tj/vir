/* grand-canonical simulation to compute the partition function
 * of all biconnected diagrams
 * simple Widom insertion with heat-bath removal */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dgrjw.h"
#include "mcutil.h"



int nmin = 1; /* the minimal order of virial coefficients */
int nmax = DG_NMAX; /* the maximal order of virial coefficients */
real mcamp = 1.5f;
double nequil = 10000000;
double nsteps = 10000000;
double ratn = 0.5; /* frequency of n-moves */
  /* this version used a simple n-move, which is on average faster than
   * configurational sampling, so we should use a larger ratn */
real rc0 = 1; /* initial rc */

double mindata = 100; /* minimal # of data points to make update */
const char *fninp = NULL; /* input file */
const char *fnout = NULL; /* output file */
const char *fnZrtmp = "Zrh.tmp"; /* temporary output file */
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

int neql = 0; /* rounds of equilibration (in which Z is updated) */
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
  argopt_add(ao, "-G",    "%d",   &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-A",    "%d",   &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-F",    "%d",   &nstfb,   "interval of computing fb");
  argopt_add(ao, "-X",    "%d",   &nedxmax, "maximal number of extra edges for computing fb");
  argopt_add(ao, "-Q",    "%lf",  &mindata, "minimal number of data points to warrant parameter update");
  argopt_add(ao, "-q",    "%d",   &nstsave, "interval of saving data");
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

  /* Metropolis move */
  if (nmvtype == 0) updrc = 0;

  /* make nstcs a multiple of nsted */
  nstcs = (nstcs + nsted - 1) / nsted * nsted;

  if (nstfb < 0) {
    if (D <= 10) {
      nstfb = (100 + nstcs - 1) / nstcs * nstcs;
    } else {
      nstfb = nstcs;
    }
  }

  if (nedxmax < 0) {
    int nedarr[] = {0, 9, 9, 9, 9, 9, 9, 11, 11, 12, 12, 12, 13, 13};
    nedxmax = (D < 14) ? nedarr[D] : 14;
  }

  argopt_dump(ao);
  argopt_close(ao);
}



/* clear data */
static void cleardata(int nmax, double *hist,
    double (*nup)[2], double (*ndown)[2], double (*ngpr)[2],
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  int i;

  for (i = 0; i <= nmax; i++) {
    hist[i] = 0;
    nup[i][0] = ndown[i][0] = 1e-20;
    nup[i][1] = ndown[i][1] = 0;
    ngpr[i][0] = 1e-20;
    ngpr[i][1] = 0;
    nedg[i][0] = 1e-20;
    nedg[i][1] = 0;
    ncsp[i][0] = 1e-20;
    ncsp[i][1] = 0;
    fbsm[i][0] = 1e-20;
    fbsm[i][1] = fbsm[i][2] = 0;
  }
}



/* update Zr */
static void updateZr(int nmin, int nmax, int updrc,
    double *Zr, double *rc,
    const double *Zr0, const double *rc0, double mindata,
    double (*nup)[2], double (*ndown)[2], double (*ngpr)[2])
{
  int i;
  double r, nZr;

  if (Zr != Zr0) for (i = 0; i <= nmax; i++) Zr[i] = Zr0[i];
  if (rc != rc0) for (i = 0; i <= nmax; i++) rc[i] = rc0[i];

  for (i = nmin + 1; i <= nmax; i++) {
    if (i == 2) continue; /* exact */
    r = getrrat(nup[i - 1][1], nup[i - 1][0], ndown[i][1], ndown[i][0],
        mindata, 0.7, 1.4);
    if (nmvtype == 1 && updrc && ngpr[i][1] > mindata) {
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
      nZr = ngpr[i][0] / ngpr[i][1];
      rc[i] *= pow(r * Zr[i] / nZr, 1. / D);
      Zr[i] = nZr;
    } else { /* update Zr only */
      Zr[i] *= r;
    }
    if (nmvtype == 0) rc[i] = pow(Zr[i], 1. / D);
  }
}



/* compute the partition function
 * assuming updateZr() has been called */
static void computeZ(int nmax, double *Z, double *B, const double *Zr,
    const double *rc, double (*fbsm)[3])
{
  int i;
  double z, fbav, fac;

  Z[0] = Z[1] = 1;
  B[0] = B[1] = 1;
  for (fac = 1, z = 1, i = 2; i <= nmax; i++) {
    /* The adjusted acceptance ratios satisfy
     *  AR(n->n-1)/AR(n-1->n) = 1 = rc*^D f* n (n - 1) / ZZ
     * where, ZZ = Z(n)/Z(n-1) is the ratio of the actual
     * partition function. */
    z *= Zr[i];
    if (nmvtype == 1) z *= pow(rc[i], D) * i * (i - 1);
    if (i <= 4) { /* reset the exact values */
      if (i <= 2) z = 1;
      else if (i == 3) {
        printf("Z3: %g to %g (ref)\n", z, Z3rat(D));
        z = Z3rat(D);
      } else if (i == 4) {
        printf("Z4: %g vs %g (ref. %s)\n", z, Z4rat(D), D > 12 ? "approx." : "exact");
        if (D < 12) z = Z4rat(D);
      }
    }
    fac *= 2./i;
    fbav = fbsm[i][1] / fbsm[i][0];
    Z[i] = z;
    /* B(i)/B(2)^(i-1) = (1 - i) 2^(i-1) / i! Z(i)/Z(2)^(i-1) <fb(i)> */
    B[i] = (1 - i) * fac * z * fbav;
  }
}



/* print arrays to screen */
static void printZr(int nmax, int nmin, const double *Zr,
    const double *rc, const double *Z,
    const double *hist, double (*nup)[2], double (*ndown)[2],
    double (*ngpr)[2], double (*nedg)[2], double (*ncsp)[2],
    double (*fbsm)[3], const double *B)
{
  int i;

  for (i = nmin; i <= nmax; i++) {
    printf("%2d %9.6f %8.6f %.6e %12.0f %.5f %.5f %7.4f %7.4f %10.8f %9.6f %+9.6f %+.7e\n",
        i, Zr[i], rc[i], Z[i], hist[i],
        nup[i-1][1] / nup[i-1][0], ndown[i][1] / ndown[i][0],
        ngpr[i][1] / ngpr[i][0],
        nedg[i][1] / nedg[i][0], ncsp[i][1] / ncsp[i][0],
        fbsm[i][2] / fbsm[i][0], fbsm[i][1] / fbsm[i][0], B[i]);
  }
}



#define mkfnZrdef(fn, fndef, d, n) if (fn == NULL) { \
  sprintf(fndef, "Zr%sD%dn%d.dat", (nmvtype == 1) ? "h" : "", d, n); \
  fn = fndef; }



/* save the partition function to file
 * we save the ratio of the successive values, such that one can modify
 * a particular value without affecting later ones */
static int saveZr(const char *fn, int nmin, int nmax,
    const double *Zr, const real *rc,
    const double *Z, const double *hist,
    double (*nup)[2], double (*ndown)[2],
    double (*ngpr)[2],
    double (*nedg)[2], double (*ncsp)[2],
    double (*fbsm)[3], const double *B)
{
  FILE *fp;
  char fndef[64];
  int i;
  double tot = 0;

  mkfnZrdef(fn, fndef, D, nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "#%s %d %d %d %d V3 %d\n", (nmvtype == 1) ? "H" : "",
      D, nmax, nmin, nmvtype, nedxmax);
  for (i = nmin; i <= nmax; i++) {
    fprintf(fp, "%3d %18.14f %20.14e %18.14f "
        "%14.0f %14.0f %14.0f %.14f %.14f ",
        i, Zr[i], Z[i], rc[i], hist[i], nup[i-1][0], ndown[i][0],
        nup[i-1][1] / nup[i-1][0], ndown[i][1] / ndown[i][0]);
    if (nmvtype == 1)
      fprintf(fp, "%14.0f %.14f ", ngpr[i][0], ngpr[i][1] / ngpr[i][0]);
    fprintf(fp, "%14.0f %14.0f %14.0f %18.14f %16.14f %16.14f %+17.14f %+20.14e\n",
        nedg[i][0], ncsp[i][0], fbsm[i][0],
        nedg[i][1] / nedg[i][0], ncsp[i][1] / ncsp[i][0],
        fbsm[i][2] / fbsm[i][0], fbsm[i][1] / fbsm[i][0], B[i]);
    tot += hist[i];
  }
  fclose(fp);
  printf("saved Zr to %s, tot %g\n", fn, tot);
  return 0;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int loadZr(const char *fn, double *Zr, real *rc, int nmax,
    int loadall, double *hist, double (*nup)[2], double (*ndown)[2],
    double (*ngpr)[2],
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  FILE *fp;
  int i = -1, d = 0, n = 0, n0 = 1, offset = 0, next = 0, ver = 0;
  char s[512], fndef[64], sver[16] = "", *p;

  mkfnZrdef(fn, fndef, D, nmax);
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
  if ( ver >= 3
    && ( sscanf(s + 2 + next, "%d", &i) != 1 || i != nmvtype ) ) {
    fprintf(stderr, "%s: mvtype %d vs %d\n%s", fn, i, nmvtype, s);
    fclose(fp);
    return -1;
  }
  if (loadall && ver <= 2) {
    fprintf(stderr, "%s: %s does not have enough data\n", fn, sver);
    loadall = 0;
  }

  offset = (ver == 0) ? 1 : 0;
  Zr[0] = Zr[1] = 1;
  if (ver == 0) n0 = 2; else if (ver <= 2) n0 = 3;
  for (i = n0; i <= nmax; i++) {
    double Z, rc1;

    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if ( 4 != sscanf(s, "%d%lf%lf%lf%n", &n, &Zr[i], &Z, &rc1, &next)
      || i != n + offset) {
      fprintf(stderr, "%s V%d: ends on line %d(n %d)\n%s", fn, ver, i, n, s);
      break;
    }
    rc[i] = (real) rc1;
    if (loadall) { /* try to get additional data */
      p = s + next;
      if ( 5 != sscanf(p, "%lf%lf%lf%lf%lf%n", &hist[i], &nup[i-1][0],
                &ndown[i][0], &nup[i-1][1], &ndown[i][1], &next) ) {
        fprintf(stderr, "%s: corrupted(I) on line %d\n%s", fn, i, s);
        break;
      }
      p += next;
      if (nmvtype == 1) {
        if ( 2 != sscanf(p, "%lf%lf%n", &ngpr[i][0], &ngpr[i][1], &next) ) {
          fprintf(stderr, "%s: corrupted(II) on line %d\n%s", fn, i, s);
          break;
        }
        p += next;
      }
      if ( 7 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf",
             &nedg[i][0], &ncsp[i][0], &fbsm[i][0],
             &nedg[i][1], &ncsp[i][1], &fbsm[i][2], &fbsm[i][1]) ) {
        fprintf(stderr, "%s: corrupted(III) on line %d\n%s", fn, i, s);
        break;
      }
      if (nup[i-1][0] <= 0) nup[i-1][0] = 1e-20;
      if (ndown[i][0] <= 0) ndown[i][0] = 1e-20;
      nup[i-1][1] *= nup[i-1][0];
      ndown[i][1] *= ndown[i][0];
      if (ngpr[i][0] <= 0) ngpr[i][0] = 1e-20;
      if (nedg[i][0] <= 0) nedg[i][0] = 1e-20;
      if (ncsp[i][0] <= 0) ncsp[i][0] = 1e-20;
      if (fbsm[i][0] <= 0) fbsm[i][0] = 1e-20;
      ngpr[i][1] *= ngpr[i][0];
      nedg[i][1] *= nedg[i][0];
      ncsp[i][1] *= ncsp[i][0];
      fbsm[i][1] *= fbsm[i][0];
      fbsm[i][2] *= fbsm[i][0];
    }
  }
  fclose(fp);
  printf("loaded Zr from %s\n", fn);
  return i - 1;
}



/* accumulate data */
static void accumdata(const dg_t *g, double t, int nstcs, int nstfb,
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  int ned = -1, n = g->n, ncs, sc, err, fb = 0;
  int nlookup = DGMAP_NMAX;
  static int degs[DG_NMAX];

  /* we will always compute fb if a lookup table is available
   * i.e., for a small n, but the n = 8 case costs 30s-70s,
   * which is a bit long for a short simulation */
  if (nsteps < 5e8) nlookup = 7;

  ned = dg_degs(g, degs);
  nedg[n][0] += 1;
  nedg[n][1] += ned;

  if (n <= nlookup || ((int) fmod(t, nstcs) == 0) ) {
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
    ncsp[n][0] += 1;
    ncsp[n][1] += ncs;

    if (n <= nlookup || (nstfb > 0 && (int) fmod(t, nstfb) == 0)) {
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
      fbsm[n][0] += 1;
      fbsm[n][1] += fb;
      fbsm[n][2] += abs(fb);
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
static void mcgc(int nmin, int nmax, double nsteps, double mcamp,
    int neql, double nequil, int nstsave)
{
  int i, j, k, n, npr, deg, Znmax, acc, it, ieql;
  dg_t *g, *ng;
  double t, cacc = 0, ctot = 0;
  rvn_t x[DG_NMAX], xi;
  real r2ij[DG_NMAX][DG_NMAX], r2i[DG_NMAX];

  /* parameters */
  double Zr[DG_NMAX + 1] = {0}; /* ratio of the partition function, measured by V */
  double Zr1[DG_NMAX + 1] = {0}; /* alternative buffer */
  real rc1[DG_NMAX + 1]; /* alternative buffer */
  double Z[DG_NMAX + 1]; /* partition function (for output) */
  double B[DG_NMAX + 1]; /* virial coefficients (for output) */
  double hist[DG_NMAX + 1]; /* histogram */
  double nedg[DG_NMAX + 1][2]; /* number of edges */
  double ncsp[DG_NMAX + 1][2]; /* no clique separator */
  double fbsm[DG_NMAX + 1][3]; /* sum of weights */
  double ndown[DG_NMAX + 1][2]; /* acceptance probabilities */
  double nup[DG_NMAX + 1][2]; /* acceptance probabilities */
  double ngpr[DG_NMAX + 1][2]; /* good pairs */
  real rc[DG_NMAX + 1];

  cleardata(nmax, hist, nup, ndown, ngpr, nedg, ncsp, fbsm);
  /* initialize the partition function */
  if (nmvtype == 1) { /* heat-bath algorithm */
    for (i = 0; i <= DG_NMAX; i++) {
      Zr[i] = 1;
      rc[i] = rc0;
    }
    Zr[2] = 0.5;
    rc[2] = 1;
  } else { /* Metropolis algorithm */
    Zr[0] = Zr[1] = Zr[2] = 1;
    Zr[3] = Z3rat(D);
    Zr[4] = Z4rat(D) / Zr[3];
    for (i = 5; i <= DG_NMAX; i++) Zr[i] = 0.66 * i;
  }
  i = loadZr(fninp, Zr, rc, nmax, restart && !bsim0,
      hist, nup, ndown, ngpr, nedg, ncsp, fbsm);
  Znmax = (i > nmax) ? i : nmax;
  for (i = nmin; i <= nmax; i++) {
    if (nmvtype == 0) rc[i] = pow(Zr[i], 1./D);
    printf("%3d: %.6f %.6f\n", i, Zr[i], rc[i]);
  }

  g = dg_open(nmax);
  ng = dg_open(nmax);
  initx(x, nmax);
  /* initial n in nmin */
  calcr2ij(r2ij, x, nmin);
  mkgraphr2ij(g, r2ij, 1, nmin);

  ieql = (neql != 0); /* rounds of equilibrations */

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    n = g->n;
    die_if (n < nmin || n > nmax, "bad n %d, t %g\n", n, t);
    if (rnd0() < ratn) { /* switching the ensemble, O(n) */
      if (rnd0() < 0.5) { /* add an vertex */
        if (n >= nmax) goto STEP_END;
        /* attach a new vertex to a random vertex */
        i = (n == 1) ? 0 : (int) (rnd0() * n);
        rvn_inc(rvn_rndball(x[n], rc[n + 1]), x[i]);
        /* test if adding x[n] leaves the diagram biconnected */
        for (deg = 0, j = 0; j < n; j++)
          if ( (r2i[j] = rvn_dist2(x[n], x[j])) < 1 )
            ++deg;
        if ( n == 1 ) deg *= 2;
        /* the extended configuration is biconnected
         * if x[n] is connected two vertices */
        nup[n][0] += 1;
        if ( deg >= 2 ) {
          g->n = n + 1;
          for (j = 0; j < n; j++) {
            r2ij[j][n] = r2ij[n][j] = r2i[j];
            if ( r2i[j] < 1)
              dg_link(g, j, n);
            else
              dg_unlink(g, j, n);
          }
          if (nmvtype == 1) { /* heat-bath algorithm */
            npr = getpair(NULL, NULL, g, r2ij, rc[n + 1]);
            if (rnd0() < 1./(Zr[n + 1] * npr)) {
              /* accept the move */
              nup[n][1] += 1;
            } else { /* recover */
              dg_remove1(g, g, n);
            }
          } else { /* Metropolis algorithm */
            nup[n][1] += 1;
          }
        }
      } else if (n > nmin) { /* remove a random vertex */
        ndown[n][0] += 1;
        if (nmvtype == 1) { /* heat-bath removal */
          if ( (npr = getpair(&i, &j, g, r2ij, rc[n])) <= 0 )
            goto STEP_END;
          ngpr[n][0] += 1;
          ngpr[n][1] += npr;
          die_if (n == 2 && npr != 2, "n %d, npr %d\n", n, npr);
          acc = (rnd0() < npr * Zr[n]);
        } else { /* Metropolis removal */
          i = randpair(n, &j);
          die_if (i == j || i >= n || j >= n, "bad i %d, j %d, n %d\n", i, j, g->n);
          /* test if * the graph is biconnected without i
           * and the pair i and j are connected */
          acc = ( ( n <= 2
                 || dg_biconnectedvs(g, mkbitsmask(n) ^ MKBIT(i)) )
              && rvn_dist2(x[i], x[j]) < rc[n] * rc[n] );
        }

        if (acc) { /* accept the removal */
          ndown[n][1] += 1;
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
    }
    else /* configuration sampling */
    {
      if (n == 1) acc = 1;
      else BCSTEPR2(acc, i, n, g, ng, x, xi, r2ij, r2i, mcamp, gaussdisp);
      cacc += acc;
      ctot += 1;
    }
STEP_END:
    hist[g->n] += 1;
    if (it % nsted == 0)
      accumdata(g, t, nstcs, nstfb, nedg, ncsp, fbsm);
    if (it % nstcom == 0) /* remove the center of mass motion */
      shiftr2ij(g, x, r2ij);
    if (ieql) { /* equilibration */
      if ((int) fmod(t + .5, nequil) == 0) {
        if (neql > 0) { /* active equilibration */
          updateZr(nmin, nmax, updrc, Zr, rc, Zr, rc, mindata, nup, ndown, ngpr);
          computeZ(nmax, Z, B, Zr, rc, fbsm);
          saveZr(fnZrtmp, nmin, Znmax, Zr, rc, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
          printZr(nmax, nmin, Zr, rc, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
          printf("equilibration stage %d/%d\n", ieql, neql);
          cleardata(nmax, hist, nup, ndown, ngpr, nedg, ncsp, fbsm);
          t = 0; it = 0; /* reset time */
          if (++ieql >= neql) /* stop equilibration */
            ieql = 0;
        } else { /* passive equilibration */
          cleardata(nmax, hist, nup, ndown, ngpr, nedg, ncsp, fbsm);
          it = 0; it = 0; /* reset time */
          printf("t %g: equilibrated, n %d\n", t, g->n);
          t = 0; it = 0; /* reset time */
          ieql = 0;
        }
      }
    } else { /* production */
      if (it % nstsave == 0 || t > nsteps - 0.5) {
        updateZr(nmin, nmax, updrc, Zr1, rc1, Zr, rc, mindata, nup, ndown, ngpr);
        computeZ(nmax, Z, B, Zr1, rc1, fbsm);
        if (restart) { /* don't write the new Zr for a restartable simulation */
          saveZr(fnout, nmin, Znmax, Zr, rc, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
          saveZr(fnZrtmp, nmin, Znmax, Zr1, rc1, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
        } else {
          saveZr(fnout, nmin, Znmax, Zr1, rc1, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
        }
        it = 0;
      }
    }
  }
  printZr(nmax, nmin, Zr1, rc1, Z, hist, nup, ndown, ngpr, nedg, ncsp, fbsm, B);
  printf("cacc %g\n", cacc/ctot);
  dg_close(g);
  dg_close(ng);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mcgc(nmin, nmax, nsteps, mcamp, neql, nequil, nstsave);
  mtsave(NULL);
  return 0;
}

