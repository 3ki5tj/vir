/* grand-canonical simulation */
#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"

int nmin = 2; /* the minimal number of particles */
int nmax = 20;
real mcamp = 1.5f;
double nequil = 100000;
double nsteps = 10000000;
double ratn = 0.1;
real rc = 1;


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
  argopt_parse(ao, argc, argv);

  /* Monte Carlo move amplitude */
  mcamp /= D;

  argopt_dump(ao);
  argopt_close(ao);
}



/* grand canonical simulation */
static void mcgc(int nmin, int nmax, real rc,
    double nsteps, double mcamp)
{
  int i, j, deg, *conn, acc;
  dg_t *g, *ng;
  double t, vol, cacc = 0, ctot = 0;
  double *Z; /* partition function */
  double *his; /* histogram */
  double *nacc; /* acceptance probabilities */
  rvn_t *x, xi;
  real *xn;

  xnew(Z, nmax + 1);
  xnew(his, nmax + 1);
  for (i = 0; i <= nmax; i++) {
    Z[i] = 1;
    his[i] = 0;
  }
  xnew(nacc, 4 * nmax + 4);
  for (i = 0; i <= nmax * 4 + 3; i++)
    nacc[i] = 1e-6;
  g = dg_open(nmax);
  ng = dg_open(nmax);
  xnew(x, nmax);
  xnew(conn, nmax);
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
        nacc[g->n*4] += 1;
        //printf("proposing..., deg %d, n %d, %g %g, rc %g, nacc %g, %g\n", deg, g->n, rvn_dist(xi, x[0]), rvn_dist(xi, x[1]), rc, nacc[g->n*4], nacc[g->n*4+1]);
        /* the extended configuration is biconnected if xi
         * is connected two vertices */
        if ( deg >= 2 && rnd0() < Z[g->n + 1] / Z[g->n] ) {
          nacc[g->n*4 + 1] += 1;
          //printf("accepting..., n %d, nacc %g, %g\n", g->n, nacc[g->n*4], nacc[g->n*4+1]);
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
        nacc[g->n*4 + 2] += 1;
        if ( rvn_dist2(x[g->n - 1], x[0]) < rc * rc
          && dg_biconnectedvs(g, (1u << (g->n - 1)) - 1)
          && rnd0() < Z[g->n - 1] / Z[g->n] )
        {
          nacc[g->n*4 + 3] += 1;
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
    his[g->n] += 1;
  }
  for (i = nmin; i <= nmax; i++)
    printf("%4d %12.0f %.8f %.8f %.8f\n", i, his[i], Z[i],
        nacc[4*i+1]/nacc[4*i], nacc[4*i+3]/nacc[4*i+2]);
  printf("cacc %g\n", cacc/ctot);
  dg_close(g);
  dg_close(ng);
  free(x);
  free(conn);
  free(Z);
  free(his);
  free(nacc);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  mcgc(nmin, nmax, rc, nsteps, mcamp);
  mtsave(NULL);
  return 0;
}

