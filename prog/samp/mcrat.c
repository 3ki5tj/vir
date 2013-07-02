/* compute the ratio of two successive virial coefficients */


/* #define D to define the dimension */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgsc.h"


typedef long long long_t;

int n = 4; /* order */
long_t nstepsmc = (long_t) 10000000ll;



/* compute the sign of the virial coefficient */
static double mcrun(int n, long_t nsteps, int plus, double *bcrat)
{
  rvn_t *x, xi, xc;
  int i, j, sc, nedg;
  long_t scsum = 0, t, acc = 0;
  double bcsum = 0, bctot = 0;
  real amp = (real) 2.0/D;
  dg_t *g, *ng, *sg;
  double scrat;

  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  dg_full(g); /* start with fully connected diagram */
  sc = dg_rhsc(g);
  nedg = dg_nedges(g);

  for (t = 0; t < nsteps; t++) {
    i = (int) (rnd0() * n);
    rvn_rnddisp(xi, x[i], amp);
    dg_copy(ng, g);
    for (j = 0; j < n; j++) {
      if (j == i) continue;
      if (rvn_dist2(xi, x[j]) >= 1)
        dg_unlink(ng, i, j);
      else
        dg_link(ng, i, j);
    }
    if ( dg_biconnected(ng) ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      sc = dg_rhsc(ng);
      nedg = dg_nedges(ng);
      acc++;
    }
    scsum += sc * (nedg % 2 ? -1 : 1);

#if 0
    printf("t %lld, i %d\n", t, i);
    for (i = 0; i < n; i++)
      for (j = i + 1; j < n; j++)
        if ((rvn_dist2(x[i], x[j]) < 1) != dg_linked(g, i, j)) {
          printf("connectivity corruptions t %lld, i %d, j %d\n",
              t, i, j);
          dg_print(g);
          dg_empty(ng);
          for (i = 0; i < n; i++)
            for (j = i + 1; j < n; j++){
              if (rvn_dist2(x[i], x[j]) < 1)
                dg_link(ng, i, j);
            }
          dg_print(ng);
          exit(1);
        }
#endif

    if (t % 5 == 0) { /* compute the acceptance rates */
      /* remove center of mass */
      rvn_zero(xc);
      for (j = 0; j < n; j++) rvn_inc(xc, x[j]);
      rvn_smul(xc, (real) 1./n);
      for (j = 0; j < n; j++) rvn_dec(x[j], xc);

      i = (int) (rnd0() * n);
      if (plus) {
        int deg = 0;
        real rad, r2, r2m = 0;
        for (j = 0; j < n; j++)
          if ((r2 = rvn_sqr(x[j])) > r2m)
            r2m = r2;
        rad = (real) (sqrt(r2m) + 1);
        /* see if add a point makes it connected */
        rvn_rndball(xi, rad);
        /* since the the new xi is connected to x[i],
         * biconnectivity == if xi is connected to another vertex */
        for (j = 0; j < n; j++)
          if (rvn_dist2(xi, x[j]) < 1)
            if (++deg >= 2) {
              bcsum += 1;
              break;
            }
        /* the denominator is the sampling volume */
        bctot += 1./(rad * rad * rad);

      } else {
        /* see if n - 1 points are biconnected */
        //if (rvn_dist2(x[0], x[i]) < 1)
        {
          dg_shrink1(sg, g, i);
          bcsum += dg_biconnected(sg);
          bctot += 1;
        }
      }
    }
  }
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);

  scrat = 1. * scsum / nsteps;
  *bcrat = 1. * bcsum / bctot;
  fprintf(stderr, "n %d, acc %.8f, sc %.8f, bc %.8f\n",
      n, 1.*acc/nsteps, scrat, *bcrat);
  return scrat;
}



int main(int argc, char **argv)
{
  double rate[2], bcrat[2], rvir;

  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nstepsmc = (long_t) atof(argv[2]);
  printf("n = %d, D = %d, nsteps = %lld\n", n, D, nstepsmc);

  rate[0] = mcrun(n, nstepsmc, 0, bcrat);
  rate[1] = mcrun(n - 1, nstepsmc, 1, bcrat + 1);
  rvir = rate[0]/rate[1] * (bcrat[1]/bcrat[0]) * (1. - n)/(2. - n)/n;
  printf("n %d, rate %g/%g, bcrat %g/%g, vir%d/vir%d %g, %g\n",
      n, rate[0], rate[1], bcrat[0], bcrat[1],
      n, n - 1, rvir, rvir * 2);

  mtsave(NULL);
  return 0;
}

