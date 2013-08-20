/* grand-canonical simulation to compute the partition function
 * of all biconnected diagrams */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"

int nmin = 3; /* the minimal number of particles */
int nmax = 20;
real mcamp = 1.5f;
double nequil = 100000;
double nsteps = 10000000;
double ratn = 0.1;
real rc = 0;
const char *fnin = NULL; /* input file */
const char *fnout = NULL; /* output file */
int nststat = 100;


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
  argopt_add(ao, "-f",   "%d",  &nststat, "interval of accumulating data");
  argopt_parse(ao, argc, argv);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  if (rc <= 0) rc = 1;
  argopt_dump(ao);
  argopt_close(ao);
}



/* save the partition function to file
 * we save the ratio of the successive values, such that one can modify
 * a particular value without affecting later ones */
static int saveZr(const char *fn, const double *Zr, int nmax,
    const double *hist, const double *nncs, const double *nedg,
    int nststat, const double *nacc)
{
  FILE *fp;
  int i;
  char fndef[64];
  double x;

  if (fn == NULL) {
    sprintf(fndef, "ZrD%d.dat", D);
    fn = fndef;
  }
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d 0\n", D, nmax);
  for (x = 1, i = 1; i < nmax; i++) {
    x *= Zr[i];
    fprintf(fp, "%2d %24.14f %24.14e %14.0f %.14f %18.14f %.14f %.14f\n",
        i, Zr[i], x, hist[i+1], nststat * nncs[i + 1] / hist[i + 1],
        nststat * nedg[i + 1] / hist[i + 1],
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
  char s[256];

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
  double Zr32, Zr43;

  /* set the exact values */
  Zr[0] = Zr[1] = 1;
  Zr[2] = Zr32 = Z3rat(D);
  Zr[3] = Zr43 = Z4rat(D)/Zr[2]; /* approximate for D > 12 */

  if (fn != NULL) {
    i0 = loadZr(fn, Zr, DG_NMAX);
    if (i0 <= 0) i0 = 4;
    /* reset the exact results */
    Zr[0] = Zr[1] = 1;
    Zr[2] = Zr32;
    if (D <= 12) Zr[3] = Zr43;
  }
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



/* grand canonical simulation */
static void mcgc(int nmin, int nmax, real rc,
    double nsteps, double mcamp)
{
  int i, j, deg, ned = 0, conn[DG_NMAX], Znmax;
  dg_t *g, *ng;
  double t, r, vol, cacc = 0, ctot = 0;
  double Zr[DG_NMAX + 1]; /* ratio of the partition function, measured by V */
  double hist[DG_NMAX + 1]; /* histogram */
  double nncs[DG_NMAX + 1]; /* number of diagrams with no clique separator */
  double nedg[DG_NMAX + 1]; /* number of edges */
  double nacc[DG_NMAX * 4 + 8]; /* acceptance probabilities */
  rvn_t x[DG_NMAX], xi;

  vol = pow(rc, D);
  Znmax = initZr(fnin, Zr, nmax); /* initialize the partition function */
  for (i = 0; i <= DG_NMAX; i++) {
    hist[i] = 1e-14;
    nncs[i] = 0;
    nedg[i] = 0;
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
        //printf("proposing..., deg %d, n %d, %g %g, rc %g, nacc %g, %g\n", deg, g->n, rvn_dist(xi, x[0]), rvn_dist(xi, x[1]), rc, nacc[g->n*4+2], nacc[g->n*4+3]);
        /* the extended configuration is biconnected if xi
         * is connected two vertices */
        if ( deg >= 2 && rnd0() < vol / Zr[g->n] ) {
          nacc[g->n*4 + 3] += 1;
          //printf("accepting..., n %d, nacc %g, %g\n", g->n, nacc[g->n*4+2], nacc[g->n*4+3]);
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
    if ((int) fmod(t, nststat) == 0) {
      ned = dg_nedges(g);
      nedg[g->n] += ned;
      if (ned == g->n) { /* ring diagram */
        nncs[g->n] += 1;
      } else {
        nncs[g->n] += (dg_cliquesep(g) == 0);
      }
    }
  }

  for (i = nmin + 1; i <= nmax; i++) {
    r = nacc[4*i - 1] / nacc[4*i - 2] / (nacc[4*i + 1] / nacc[4*i]);
    if (i >= 5 || (i == 4 && D > 12))
      Zr[i - 1] *= r;
  }

  for (i = 1; i <= nmax; i++) {
    printf("%4d %16.8f %12.0f %.8f %12.8f %.8f %.8f\n", i, Zr[i], hist[i],
        nncs[i] / hist[i] * nststat, nedg[i] / hist[i] * nststat,
        nacc[4*i + 3] / nacc[4*i + 2], nacc[4*i + 5] / nacc[4*i + 4]);
  }
  saveZr(fnout, Zr, Znmax, hist, nncs, nedg, nststat, nacc);
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

