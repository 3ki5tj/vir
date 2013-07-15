/* compute the ratio of two successive virial coefficients
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcrat0.c -lm
 * define MVALL to use the all-particle move */


#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"



int n = 4; /* order */
double nstepsmc = 10000000.0;



/* compute a Ree-Hoover diagram from the coordinates */
static void mkgraph(dg_t *g, rvn_t *x)
{
  int i, j, n = g->n;

  dg_empty(g);
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++)
      if (rvn_dist2(x[i], x[j]) < 1)
        dg_link(g, i, j);
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nsteps, real amp, int plus, double *nrat)
{
  rvn_t *x, xi, xc;
#ifdef MVALL /* moving all particles */
  rvn_t *y;
#endif
  int i, j, fb;
  double fbav, fbsum = 0, t, acc = 0, nsum = 0, ntot = 0;
  dg_t *g, *ng, *sg;

  xnew(x, n + 1);
  amp *= (real) 1.0/D;
#ifdef MVALL
  xnew(y, n + 1);
  amp *= 0.2/sqrt(D);
#endif
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb(g);

  for (t = 1; t <= nsteps; t++) {
#ifdef MVALL
    for (i = 0; i < n; i++)
      rvn_rnddisp(y[i], x[i], amp);
    mkgraph(ng, y);
#else
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
#endif
    if ( dg_biconnected(ng) ) { /* accept the move */
#ifdef MVALL
      for (i = 0; i < n; i++)
        rvn_copy(x[i], y[i]);
#else
      rvn_copy(x[i], xi);
#endif
      dg_copy(g, ng);
      fb = dg_hsfb(g);
      acc++;
    }
    fbsum += fb;

    if (fmod(t, 10) < 0.1) { /* compute the acceptance rates */
      /* remove center of mass */
      rvn_zero(xc);
      for (j = 0; j < n; j++) rvn_inc(xc, x[j]);
      rvn_smul(xc, (real) 1./n);
      for (j = 0; j < n; j++) rvn_dec(x[j], xc);

      i = (int) (rnd0() * n);
      if (plus) {
        /* compute the probability of adding a vertex
         * that leaves the graph biconnected */
        int deg = 0;
        real rad, r2, r2m = 0, r2n = 0, vol;
        
        /* naively, in a trial, the vertex should be added uniformly to
         * a very large volume Vmax, and the probability is computed as
         *  r = (V_bi / Vmax) * (Vmax / Vunit)             (1)
         * where Vmax/Vunit is the normalization to a uniform sphere
         * and B2 = Vunit / 2; (V_bi / Vmax) is the acceptance ratio.
         * We can however the shrink the sampling volume to a sphere
         * with the radius being the second largest distance from any
         * vertex to the center of mass (because there should be at
         * least two vertices connected to the new vertex to make the
         * graph biconnected). Thus
         *  r = (V_bi / V_samp) * (Vsamp / Vunit)         (2)
         * and (V_bi / V_samp) is the new and enlarged acceptance
         * probability. */
        for (j = 0; j < n; j++)
          if ((r2 = rvn_sqr(x[j])) > r2m) {
            r2n = r2m; /* second largest volume */
            r2m = r2;
          } else if (r2 > r2n) {
            r2n = r2;
          }
        if (r2n <= 0) r2n = r2m;

        /* we only need to sample the second largest sphere */
        rad = (real) (sqrt(r2n) + 1);
        /* compute the relative volume to the unit sphere */
        for (vol = 1, j = 0; j < D; j++) vol *= rad;
        /* see if adding a vertex leaves the graph biconnected */
        rvn_rndball(xi, rad);
        /* biconnectivity means xi is connected two vertices */
        for (j = 0; j < n; j++)
          if (rvn_dist2(xi, x[j]) < 1)
            if (++deg >= 2) {
              nsum += vol; /* see Eq. (2) above */
              break;
            }
        ntot += 1;
      } else {
        /* compute the probability of biconnectivity after removing a random vertex */
        dg_shrink1(sg, g, i);
        nsum += dg_biconnected(sg);
        ntot += 1;
      }

      if (fmod(t, 100000000) < 0.1)
        fprintf(stderr, "n %d, %g steps, acc %.8f, fb %.8f, nc %.8f\n",
        n, t, acc/t, fbsum/t, nsum/ntot);
    }
  }
  free(x);
#ifdef MVALL
  free(y);
#endif
  dg_close(g);
  dg_close(ng);
  dg_close(sg);

  fbav = 1. * fbsum / nsteps;
  *nrat = 1. * nsum / ntot;
  fprintf(stderr, "n %d, acc %.8f, fb %.8f, nc %.8f\n",
      n, 1.*acc/nsteps, fbav, *nrat);
  return fbav;
}



int main(int argc, char **argv)
{
  double fbav[2], nrat[2], rvir, amp = 1.0;

  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nstepsmc = atof(argv[2]);
  if (argc > 3) amp = atof(argv[3]);
  printf("n %d, D %d, nsteps %g, amp %g, code %d-bit\n",
      n, D, nstepsmc, amp, (int) sizeof(code_t) * 8);

  fbav[0] = mcrun(n, nstepsmc, (real) amp, 0, nrat);
  fbav[1] = mcrun(n - 1, nstepsmc, (real) amp, 1, nrat + 1);
  rvir = fbav[0]/fbav[1] * (nrat[1]/nrat[0]) * (1. - n)/(2. - n)/n;
  printf("n %d, fbav %g/%g, nrat %g/%g = %g, vir%d/vir%d %g, %g\n",
      n, fbav[0], fbav[1], nrat[0], nrat[1], nrat[0]/nrat[1],
      n, n - 1, rvir, rvir * 2);

  mtsave(NULL);
  return 0;
}

