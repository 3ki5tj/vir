/* grand-canonical simulation to compute the partition function
 * of all biconnected diagrams */
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

double mindata = 100; /* minimal # of data points to make update */
const char *fninp = NULL; /* input file */
const char *fnout = NULL; /* output file */
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



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "-z",    "%d",   &nmin,    "minimal order of virial coefficients");
  argopt_add(ao, "-n",    "%d",   &nmax,    "maximal order of virial coefficients");
  argopt_add(ao, "-1",    "%lf",  &nsteps,  "number of simulation steps");
  argopt_add(ao, "-a",    "%r",   &mcamp,   "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-r",    "%lf",  &ratn,    "rate of particle moves");
  argopt_add(ao, "-i",    NULL,   &fninp,    "input file");
  argopt_add(ao, "-o",    NULL,   &fnout,   "output file");
  argopt_add(ao, "-G",    "%d",   &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-A",    "%d",   &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-F",    "%d",   &nstfb,   "interval of computing fb");
  argopt_add(ao, "-X",    "%d",   &nedxmax, "maximal number of extra edges for computing fb");
  argopt_add(ao, "-Q",    "%lf",  &mindata, "minimal number of data points to warrant parameter update");
  argopt_add(ao, "-q",    "%d",   &nstsave, "interval of saving data");
  argopt_add(ao, "-E",    "%d",   &neql,    "number of equilibration stages, 0: disable, > 0: update Zr & sr, < 0: 1 round no update");
  argopt_add(ao, "-0",    "%lf",  &nequil,  "number of equilibration steps per round");
  argopt_add(ao, "-M",    "%d",   &nstcom,  "interval of centering the structure");
  argopt_add(ao, "-R",    "%b",   &restart, "run a restartable simulation");
  argopt_add(ao, "-B",    "%b",   &bsim0,   "start a restartable simulation");

  argopt_parse(ao, argc, argv);

  /* Monte Carlo move amplitude */
  mcamp /= D;

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
    double (*nup)[2], double (*ndown)[2],
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  int i;

  for (i = 0; i <= nmax; i++) {
    hist[i] = 0;
    nup[i][0] = ndown[i][0] = 1e-20;
    nup[i][1] = ndown[i][1] = 0;
  }
  for (i = 0; i <= nmax; i++) {
    nedg[i][0] = 1e-20;
    nedg[i][1] = 0;
    ncsp[i][0] = 1e-20;
    ncsp[i][1] = 0;
    fbsm[i][0] = 1e-20;
    fbsm[i][1] = fbsm[i][2] = 0;
  }
}



/* update Zr */
static void updateZr(int nmin, int nmax, double *Zr,
    double (*nup)[2], double (*ndown)[2], double mindata)
{
  int i;
  double r;

  for (i = nmin + 1; i <= nmax; i++) {
    r = (nup[i-1][1] / nup[i-1][0]) / (ndown[i][1] / ndown[i][0]);
    /* make sure enough data points and not to change the exact data */
    if ( (i >= 5 || (i == 4 && D > 12)) &&
        nup[i-1][1] >= mindata && ndown[i][1] >= mindata)
      Zr[i] *= r;
  }
}



/* compute the partition function
 * assuming updateZr() has been called */
static void computeZ(int nmax, double *Z, double *B, const double *Zr,
    double (*fbsm)[3])
{
  int i;
  double z, fbav, fac;

  Z[0] = Z[1] = 1;
  B[0] = B[1] = 1;
  for (fac = 1, z = 1, i = 2; i <= nmax; i++) {
    z *= Zr[i];
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
    /* B(i)/B(2)^(i-1) = (1-i) 2^(i-1) / i! Z(i)/Z(2)^(i-1) <fb(i)> */
    B[i] = (1 - i) * fac * z * fbav;
  }
}



/* print arrays to screen */
static void printZr(int nmax, const double *Zr, const double *Z,
    const double *hist, double (*nup)[2], double (*ndown)[2],
    double (*nedg)[2], double (*ncsp)[2],
    double (*fbsm)[3], const double *B)
{
  int i;

  for (i = nmin; i <= nmax; i++) {
    printf("%2d %9.6f %.6e %12.0f %.5f %.5f %7.4f %10.8f %9.6f %+9.6f %+.7e\n",
        i, Zr[i], Z[i], hist[i],
        nup[i-1][1] / nup[i-1][0], ndown[i][1] / ndown[i][0],
        nedg[i][1] / nedg[i][0], ncsp[i][1] / ncsp[i][0],
        fbsm[i][2] / fbsm[i][0], fbsm[i][1] / fbsm[i][0], B[i]);
  }
}



#define mkfnZrdef(fn, fndef, d, n) \
  if (fn == NULL) { sprintf(fndef, "ZrD%dn%d.dat", d, n); fn = fndef; }



/* save the partition function to file
 * we save the ratio of the successive values, such that one can modify
 * a particular value without affecting later ones */
static int saveZr(const char *fn, int nmin, int nmax,
    const double *Zr, const double *Z, const double *hist,
    double (*nup)[2], double (*ndown)[2],
    double (*nedg)[2], double (*ncsp)[2],
    double (*fbsm)[3], const double *B)
{
  FILE *fp;
  int i;
  char fndef[64];
  double tot = 0;

  mkfnZrdef(fn, fndef, D, nmax);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d V3 %d 0 %d\n", D, nmax, nmin, nedxmax);
  for (i = nmin; i <= nmax; i++) {
    fprintf(fp, "%3d %18.14f %20.14e %14.0f %14.0f %14.0f %.14f %.14f "
        "%14.0f %14.0f %14.0f %18.14f %16.14f %16.14f %+17.14f %+20.14e\n",
        i, Zr[i], Z[i], hist[i], nup[i-1][0], ndown[i][0],
        nup[i-1][1] / nup[i-1][0], ndown[i][1] / ndown[i][0],
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
static int loadZr(const char *fn, double *Zr, int nmax,
    int loadall, double *hist, double (*nup)[2], double (*ndown)[2],
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  FILE *fp;
  int i = -1, d = 0, n = 0, n0 = 3, offset = 0, next = 0, ver = 0;
  char s[512], fndef[64], sver[16] = "", *p;

  mkfnZrdef(fn, fndef, D, nmax);
  xfopen(fp, fn, "r", return -1);

  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if (s[0] == '#') {
    if (sscanf(s + 1, "%d%d%8s%d", &d, &n, sver, &n0) != 4
     || d != D) {
      fprintf(stderr, "%s dimension %d vs. D %d\n%s", fn, d, D, s);
      fclose(fp);
      return -1;
    }
    ver = atoi(sver[0] == 'V' ? sver + 1 : sver); /* strip V */
  } else { /* no info. line, give back the first line */
    rewind(fp);
  }
  if (loadall && ver <= 1) {
    fprintf(stderr, "%s version %s does not have enough data\n", fn, sver);
    loadall = 0;
  }

  /* offset is 1 in version 0 */
  offset = (ver == 0) ? 1 : 0;
  Zr[0] = Zr[1] = 1;
  if (ver == 0) n0 = 2; else if (ver <= 2) n0 = 3;
  for (i = n0; i <= nmax; i++) {
    double Z;

    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if (3 != sscanf(s, "%d%lf%lf%n", &n, &Zr[i], &Z, &next)
        || i != n + offset) {
      fprintf(stderr, "%s ends on line %d, version %d\n%s", fn, i, ver, s);
      break;
    }
    if (loadall) { /* try to get additional data */
      p = s + next;
      if (12 != sscanf(p, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
            &hist[i], &nup[i-1][0], &ndown[i][0],
            &nup[i-1][1], &ndown[i][1],
            &nedg[i][0], &ncsp[i][0], &fbsm[i][0],
            &nedg[i][1], &ncsp[i][1], &fbsm[i][2], &fbsm[i][1]) ) {
        fprintf(stderr, "%s corrupted on line %d\n%s", fn, i, s);
        break;
      }
      if (nup[i-1][0] <= 0) nup[i-1][0] = 1e-20;
      if (ndown[i][0] <= 0) ndown[i][0] = 1e-20;
      nup[i-1][1] *= nup[i-1][0];
      ndown[i][1] *= ndown[i][0];
      if (nedg[i][0] <= 0) nedg[i][0] = 1e-20;
      if (ncsp[i][0] <= 0) ncsp[i][0] = 1e-20;
      if (fbsm[i][0] <= 0) fbsm[i][0] = 1e-20;
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



/* initialize the partition function */
static int initZr(double *Zr, int nmax, int i0)
{
  int i, imax;

  /* reset the exact results */
  Zr[0] = Zr[1] = Zr[2] = 1;
  Zr[3] = Z3rat(D);
  if (D <= 12 || i0 < 4) Zr[4] = Z4rat(D)/Zr[3];
  if (i0 < 4) i0 = 4;
  /* to avoid bad data */
  for (i = 0; i <= DG_NMAX; i++)
    if (Zr[i] <= 0) Zr[i] = 1;

  if (nmax > i0) { /* guess the unknown */
    fprintf(stderr, "guessing Zr for n = %d to %d\n", i0 + 1, nmax);
    for (i = i0 + 1; i <= nmax; i++) Zr[i] = 0.66 * i;
    imax = nmax;
  } else {
    imax = i0;
  }
  fprintf(stderr, "initial Zr (i0 %d, imax %d, nmax %d):\n", i0, imax, nmax);
  for (i = 1; i <= imax; i++)
    fprintf(stderr, "%4d %16.8f\n", i, Zr[i]);
  return imax;
}



/* accumulate data */
static void accumdata(const dg_t *g, double t, int nstcs, int nstfb,
    double (*nedg)[2], double (*ncsp)[2], double (*fbsm)[3])
{
  int ned = -1, n = g->n, ncs, err;
  double fb = 0, sc;
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
        ncs = (fabs(sc) > 1e-3);
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
      fbsm[n][2] += fabs(fb);
    }
  }
}



/* grand canonical simulation, main function */
static void mcgc(int nmin, int nmax, double nsteps, double mcamp,
    int neql, double nequil, int nstsave)
{
  int i, j, n, deg, conn[DG_NMAX], Znmax, acc, it, ieql;
  dg_t *g, *ng;
  double t, cacc = 0, ctot = 0;
  rvn_t x[DG_NMAX], xi;

  /* parameters */
  double Zr[DG_NMAX + 1] = {0}; /* ratio of the partition function, measured by V */
  double Zr1[DG_NMAX + 1] = {0}; /* alternative buffer */
  double Z[DG_NMAX + 1]; /* partition function (for output) */
  double B[DG_NMAX + 1]; /* virial coefficients (for output) */
  double hist[DG_NMAX + 1]; /* histogram */
  double nedg[DG_NMAX + 1][2]; /* number of edges */
  double ncsp[DG_NMAX + 1][2]; /* no clique separator */
  double fbsm[DG_NMAX + 1][3]; /* sum of weights */
  double ndown[DG_NMAX + 1][2]; /* acceptance probabilities */
  double nup[DG_NMAX + 1][2]; /* acceptance probabilities */
  real rc[DG_NMAX + 1], rc2[DG_NMAX + 1], vol[DG_NMAX + 1];

  cleardata(nmax, hist, nup, ndown, nedg, ncsp, fbsm);
  i = loadZr(fninp, Zr, nmax, restart && !bsim0,
      hist, nup, ndown, nedg, ncsp,fbsm);
  Znmax = initZr(Zr, nmax, i); /* initialize the partition function */
  for (i = 1; i <= nmax; i++) {
    vol[i] = Zr[i]; /* setting vol == Zr avoid the additional
                       acceptance step */
    rc[i] = pow(vol[i], 1./D);
    rc2[i] = rc[i] * rc[i];
  }
  g = dg_open(nmax);
  ng = dg_open(nmax);
  initx(x, nmax);
  mkgraph(g, x, nmin);

  ieql = (neql > 0); /* rounds of equilibrations */

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
        for (deg = 0, j = 0; j < n; j++) {
          if ( rvn_dist2(x[n], x[j]) < 1 ) {
            conn[j] = 1;
            ++deg;
          } else {
            conn[j] = 0;
          }
        }
        if ( n == 1 ) deg *= 2;
        /* the extended configuration is biconnected
         * if x[n] is connected two vertices */
        nup[n][0] += 1;
        if ( deg >= 2 ) {
          nup[n][1] += 1;
          g->n  = n + 1;
          for (j = 0; j < n; j++) {
            if ( conn[j] )
              dg_link(g, j, n);
            else
              dg_unlink(g, j, n);
          }
        }
      } else if (n > nmin) { /* remove a random vertex */
        ndown[n][0] += 1;
        i = randpair(n, &j); /* random root at j */
        die_if (i == j || i >= n || j >= n, "bad i %d, j %d, n %d\n", i, j, g->n);
        /* test if * the graph is biconnected without i
         * and the pair i and j are connected */
        if ( ( n == 2 || dg_biconnectedvs(g, mkbitsmask(n) ^ MKBIT(i)) )
          && rvn_dist2(x[i], x[j]) < rc2[n] ) {
          ndown[n][1] += 1;
          //for (j = i; j < n - 1; j++) rvn_copy(x[j], x[j + 1]);
          dg_remove1(g, g, i); /* g->n is decreased by 1 here */
          if (i < n - 1)
            memmove(x[i], x[i + 1], (n - 1 - i) * sizeof(rvn_t));
        }
      }
    }
    else /* configuration sampling */
    {
      if ( n > 1 ) BCSTEP(acc, i, n, g, ng, x, xi, mcamp, gaussdisp);
      ctot += 1;
      cacc += acc;
    }
STEP_END:
    hist[g->n] += 1;
    if (it % nsted == 0)
      accumdata(g, t, nstcs, nstfb, nedg, ncsp, fbsm);
    if (it % nstcom == 0) /* remove the center of mass motion */
      rvn_rmcom(x, g->n);
    if (ieql) { /* equilibration */
      if ((int) fmod(t + .5, nequil) == 0) {
        updateZr(nmin, nmax, Zr, nup, ndown, mindata);
        computeZ(nmax, Z, B, Zr, fbsm);
        saveZr("Zr.tmp", nmin, Znmax, Zr, Z, hist, nup, ndown, nedg, ncsp, fbsm, B);
        printZr(nmax, Zr, Z, hist, nup, ndown, nedg, ncsp, fbsm, B);
        printf("equilibration stage %d/%d\n", ieql, neql);
        cleardata(nmax, hist, nup, ndown, nedg, ncsp, fbsm);
        t = 0; it = 0; /* reset time */
        if (++ieql >= neql) /* stop equilibration */
          ieql = 0;
      }
    } else { /* production */
      if (it % nstsave == 0 || t > nsteps - 0.5) {
        memcpy(Zr1, Zr, sizeof(Zr[0]) * (nmax + 1));
        updateZr(nmin, nmax, Zr1, nup, ndown, mindata);
        computeZ(nmax, Z, B, Zr1, fbsm);
        if (restart) { /* don't write the new Zr for a restartable simulation */
          saveZr(fnout, nmin, Znmax, Zr, Z, hist, nup, ndown, nedg, ncsp, fbsm, B);
        } else {
          saveZr(fnout, nmin, Znmax, Zr1, Z, hist, nup, ndown, nedg, ncsp, fbsm, B);
        }
        it = 0;
      }
    }
  }
  printZr(nmax, Zr1, Z, hist, nup, ndown, nedg, ncsp, fbsm, B);
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

