/* grand-canonical simulation to compute the partition function
 * of all biconnected diagrams */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dgrjw.h"
#include "mcutil.h"



int nmin = 3; /* the minimal number of particles */
int nmax = DG_NMAX - 1;
real mcamp = 1.5f;
double nequil = 100000;
double nsteps = 10000000;
double ratn = 0.1;
real rc = 0;
const char *fnin = NULL; /* input file */
const char *fnout = NULL; /* output file */
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


/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "--n0", "%d",  &nmin,    "minimal order n");
  argopt_add(ao, "-n",   "%d",  &nmax,    "maximal order n");
  argopt_add(ao, "-0",   "%lf", &nequil,  "number of equilibration steps");
  argopt_add(ao, "-1",   "%lf", &nsteps,  "number of simulation steps");
  argopt_add(ao, "-a",   "%r",  &mcamp,   "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-r",   "%lf", &ratn,    "rate of particle moves");
  argopt_add(ao, "-c",   "%r",  &rc,      "radius of particle insersion");
  argopt_add(ao, "-i",   NULL,  &fnin,    "input file");
  argopt_add(ao, "-o",   NULL,  &fnout,   "output file");
  argopt_add(ao, "-e",   "%d",  &nsted,   "interval of computing the # of edges");
  argopt_add(ao, "-s",   "%d",  &nstcs,   "interval of computing clique separator");
  argopt_add(ao, "-f",   "%d",  &nstfb,   "interval of computing fb");
  argopt_add(ao, "-x",   "%d",  &nedxmax, "maximal number of extra edges for computing fb");
  argopt_parse(ao, argc, argv);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  if (rc <= 0) rc = 1;

  if (nstfb < 0) {
    if (D <= 10) {
      nstfb = (100 + nstcs - 1) / nstcs * nstcs;
    } else {
      nstfb = nstcs;
    }
  }

  if (nedxmax < 0) {
    if (D < 8) nedxmax = 10;
    else if (D < 15) nedxmax = 12;
    else nedxmax = 14; /* a generous upper bound */
  }
  argopt_dump(ao);
  argopt_close(ao);
}



#define mkfnZrdef(fn, fndef, d) \
  if (fn == NULL) { sprintf(fndef, "ZrD%d.dat", D); fn = fndef; }



/* save the partition function to file
 * we save the ratio of the successive values, such that one can modify
 * a particular value without affecting later ones */
static int saveZr(const char *fn, const double *Zr, int nmax,
    const double *hist, const double *ncsp, const double *fbsm,
    const double *nedg, const double *nacc)
{
  FILE *fp;
  int i;
  char fndef[64];
  double x, fbav, fact;

  mkfnZrdef(fn, fndef, D);
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d 0 %d\n", D, nmax, nedxmax);
  for (fact = 1, x = 1, i = 1; i < nmax; i++) {
    x *= Zr[i];
    fact *= 2./(i + 1);
    fbav = fbsm[3*i + 4] / fbsm[3*i + 3];
    fprintf(fp, "%3d %18.14f %20.14e %14.0f "
        "%16.14f %+17.14f %+17.14f %16.14f %+20.14f "
        "%18.14f %.14f %.14f\n",
        i, Zr[i], x, hist[i + 1],
        ncsp[3*i + 4] / ncsp[3*i + 3], ncsp[3*i + 5] / ncsp[3*i + 3],
        fbav, fbsm[3*i + 5] / fbsm[3*i + 3],
        (-i)*fact * x * fbav, /* B(i+1) = 2^i (-i) / (i+1)! Z(i+1) <fb(i+1)> */ 
        nedg[2*i + 3] / nedg[2*i + 2],
        nacc[4*i + 3] / nacc[4*i + 2], nacc[4*i + 5] / nacc[4*i + 4]);
  }
  fclose(fp);
  return 0;
}



/* load the partition function from file
 * return the maximal `n' loaded */
static int loadZr(const char *fn, double *Zr, int nmax)
{
  FILE *fp;
  int i = -1, n;
  char s[512], fndef[64];

  mkfnZrdef(fn, fndef, D);
  xfopen(fp, fn, "r", return -1);

  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if (s[0] == '#') {
    if (sscanf(s + 1, "%d", &n) != 1 || n != D) {
      fprintf(stderr, "%s dimension %d vs. D %d\n%s", fn, n, D, s);
      fclose(fp);
      return -1;
    }
  } else { /* no info. line, give back the first line */
    rewind(fp);
  }

  Zr[1] = 1;
  for (i = 1; i < nmax; i++) {
    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if (2 != sscanf(s, "%d%lf", &n, &Zr[i]) || i != n) {
      fprintf(stderr, "%s ends on line %d\n%s", fn, i, s);
      break;
    }
  }
  fclose(fp);
  return i;
}



/* initialize the partition function */
static int initZr(const char *fn, double *Zr, int nmax)
{
  int i, i0 = 4, imax;

  i0 = loadZr(fn, Zr, DG_NMAX);
  /* reset the exact results */
  Zr[0] = Zr[1] = 1;
  Zr[2] = Z3rat(D);
  if (D <= 12 || i0 < 4) Zr[3] = Z4rat(D)/Zr[2];
  if (i0 < 4) i0 = 4;

  if (nmax > i0) { /* guess the unknown */
    fprintf(stderr, "guessing Zr for n = %d to %d\n", i0 + 1, nmax);
    for (i = i0; i < nmax; i++) Zr[i] = Zr[i0 - 1];
    imax = nmax;
  } else {
    imax = i0;
  }
  fprintf(stderr, "initial Zr (i0 %d, imax %d, nmax %d):\n", i0, imax, nmax);
  for (i = 1; i < imax; i++)
    fprintf(stderr, "%4d %16.8f\n", i, Zr[i]);
  return imax;
}



/* accumulate data */
static void accumdata(const dg_t *g, double t, int nsted, int nstcs, int nstfb,
    double *nedg, double *ncsp, double *fbsm)
{
  int ned = -1, n = g->n, ncs, sc, err, fb = 0;
  int nlookup = DGMAP_NMAX;
  static int degs[DG_NMAX];

  /* we will always compute fb if a lookup table is available
   * i.e., for a small n, but the n = 8 case costs 30s-70s,
   * which is a bit long for a short simulation */
  if (nsteps < 5e8) nlookup = 7;

  if (nsted > 0 && (int) fmod(t, nsted) == 0) {
    ned = dg_degs(g, degs);
    nedg[2*n] += 1;
    nedg[2*n + 1] += ned;

    if (n <= nlookup || (nstcs > 0 && (int) fmod(t, nstcs) == 0) ) {
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
          fb = sc * (1 - (ned % 2) * 2);
        } else { /* no clique separator */
          ncs = 1;
        }
      }
      ncsp[3*n] += 1;
      ncsp[3*n + 1] += ncs;
      /* this is assuming that SC = 1 for all diagrams
       * without clique separators */
      ncsp[3*n + 2] += ncs * (1 - (ned % 2) * 2);

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
        fbsm[3*n] += 1;
        fbsm[3*n + 1] += fb;
        fbsm[3*n + 2] += abs(fb);
      }
    }
  }
}



/* grand canonical simulation, main function */
static void mcgc(int nmin, int nmax, real rc,
    double nsteps, double mcamp)
{
  int i, j, deg, conn[DG_NMAX], Znmax;
  dg_t *g, *ng;
  double t, r, vol, cacc = 0, ctot = 0;
  double Zr[DG_NMAX + 1]; /* ratio of the partition function, measured by V */
  double hist[DG_NMAX + 1]; /* histogram */
  double nedg[DG_NMAX * 2 + 2]; /* number of edges */
  double ncsp[DG_NMAX * 3 + 3]; /* no clique separator */
  double fbsm[DG_NMAX * 3 + 3]; /* sum of weights */
  double nacc[DG_NMAX * 4 + 8]; /* acceptance probabilities */
  rvn_t x[DG_NMAX], xi;

  vol = pow(rc, D);
  Znmax = initZr(fnin, Zr, nmax); /* initialize the partition function */
  for (i = 0; i <= DG_NMAX; i++) {
    hist[i] = 1e-14;
    ncsp[3*i] = 1e-14;
    ncsp[3*i + 1] = ncsp[3*i + 2] = 0;
    fbsm[3*i] = 1e-14;
    fbsm[3*i + 1] = fbsm[3*i + 2] = 0;
    nedg[2*i] = 1e-14;
    nedg[2*i + 1] = 0;
  }
  for (i = 0; i < DG_NMAX * 4 + 8; i++)
    nacc[i] = 1e-6;
  g = dg_open(nmax);
  ng = dg_open(nmax);
  for (i = 0; i < nmax; i++) rvn_zero(x[i]);
  g->n = nmin;
  mkgraph(g, x);

  for (t = 1; t <= nsteps; t += 1) {
    die_if (g->n < nmin || g->n > nmax, "bad n %d, t %g\n", g->n, t);
    if (rnd0() < ratn) { /* switching the ensemble */
      if (rnd0() < 0.5) { /* add an vertex */
        if (g->n >= nmax) goto STEP_END;
        /* attach a new vertex to the first vertex */
        rvn_rndball(xi, rc);
        rvn_inc(xi, x[0]);
        /* test if adding xi leaves the diagram biconnected */
        for (deg = 0, j = 0; j < g->n; j++) {
          if ( rvn_dist2(xi, x[j]) < 1 ) {
            conn[j] = 1;
            ++deg;
          } else {
            conn[j] = 0;
          }
        }
        nacc[g->n*4 + 2] += 1;
        /* the extended configuration is biconnected if xi
         * is connected two vertices */
        if ( deg >= 2 && rnd0() < vol / Zr[g->n] ) {
          nacc[g->n*4 + 3] += 1;
          rvn_copy(x[g->n], xi);
          g->n += 1;
          for (j = 0; j < g->n - 1; j++) {
            if ( conn[j] )
              dg_link(g, j, g->n - 1);
            else
              dg_unlink(g, j, g->n - 1);
          }
        }
      } else if (g->n > nmin) { /* remove the last vertex */
        nacc[g->n*4] += 1;
        if ( rvn_dist2(x[g->n - 1], x[0]) < rc * rc
          && dg_biconnectedvs(g, (1u << (g->n - 1)) - 1)
          && rnd0() < Zr[g->n - 1] / vol)
        {
          nacc[g->n*4 + 1] += 1;
          g->n--;
          for (j = 0; j < g->n; j++)
            g->c[j] &= (1u << g->n) - 1;
        }
      }
    }
    else /* configuration sampling */
    {
      i = (int) (rnd0() * g->n);
      rvn_rnddisp(xi, x[i], mcamp);
      ng->n = g->n;
      dg_copy(ng, g);
      for (j = 0; j < g->n; j++) {
        if (j == i) continue;
        if (rvn_dist2(xi, x[j]) < 1)
          dg_link(ng, i, j);
        else
          dg_unlink(ng, i, j);
      }

      ctot += 1;
      if ( dg_biconnected(ng) ) {
        rvn_copy(x[i], xi);
        dg_copy(g, ng);
        cacc += 1;
      }
    }
STEP_END:
    hist[g->n] += 1;
    accumdata(g, t, nsted, nstcs, nstfb, nedg, ncsp, fbsm);
  }

  for (i = nmin + 1; i <= nmax; i++) {
    r = nacc[4*i - 1] / nacc[4*i - 2] / (nacc[4*i + 1] / nacc[4*i]);
    /* make sure enough data and not to erase the exact data */
    if ( (i >= 5 || (i == 4 && D > 12)) &&
        nacc[4*i - 1] >= 100 && nacc[4*i + 1] >= 100)
      Zr[i - 1] *= r;
  }

  for (i = 1; i <= nmax; i++) {
    printf("%2d %9.6f%13.0f %10.8f %+11.8f %+11.8f %11.8f %7.4f %.5f %.5f\n",
        i, Zr[i], hist[i],
        ncsp[3*i + 1] / ncsp[3*i], ncsp[3*i + 2] / ncsp[3*i],
        fbsm[3*i + 1] / fbsm[3*i], fbsm[3*i + 2] / fbsm[3*i],
        nedg[2*i + 1] / nedg[2*i],
        nacc[4*i + 3] / nacc[4*i + 2], nacc[4*i + 5] / nacc[4*i + 4]);
  }
  saveZr(fnout, Zr, Znmax, hist, ncsp, fbsm, nedg, nacc);
  printf("cacc %g\n", cacc/ctot);
  dg_close(g);
  dg_close(ng);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mcgc(nmin, nmax, rc, nsteps, mcamp);
  mtsave(NULL);
  return 0;
}

