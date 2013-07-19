/* compute the ratio of two successive virial coefficients
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcrat.c -lm */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"



int n = 8; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000, nsteps2 = -1;
real mcamp[2] = {1.5f, 0.9f};
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 100; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */
double rat10 = 0.03;
int lookup = -1; /* if to use the lookup table */


/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-2", "%lf", &nsteps2, "number of simulation steps, for n - 1");
  argopt_add(ao, "-a", "%r", &mcamp[0], "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-A", "%r", &mcamp[1], "MC amplitude for nonzero biconnected diagrams");
  argopt_add(ao, "-L", "%d", &lookup, "lookup table (-1: automatic)");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-z", "%d", &nstnmv, "interval of the n-move");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-Z", "%lf", &rat10, "estimated ratio of nonzero biconnected diagrams");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
  
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);
  if (nstfb <= 0) {
    if (lookup) {
      nstfb = 1;
    } else {
      int i;
      for (nstfb = 10, i = 8; i < n; i++) nstfb *= 3;
    }
  }
  mcamp[0] /= D;
  mcamp[1] /= D;
  /* the smaller system converges faster, so we spend less time on it */
  if (nsteps2 <= 0) nsteps2 = nsteps / 4;
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps,
    double amp[2], int nstfb, int lookup,
    int nstnmv, double *nrat)
{
  rvn_t *x, xi;
  int i, j, it, fb = 0, nfb = 0, nz, nnz, nbc, eql;
  double fbav, rathis, t;
  dg_t *g, *ng, *sg;
  double wtot[2] = {1e-6, 1e-6}, fbsum[2] = {0, 0}, fbnz[2] = {1e-6, 1e-6};
  double hist[2] = {1e-6, 1e-6}, cacc[2] = {0, 0};
  int sys, acc1;
    
  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = lookup ? dg_hsfb(g) : dg_hsfbrjw(g);
  nz = 1; /* if (lookup), nz means fb != 0,
             otherwise nz means no clique separator */
  sys = 0;
  eql = 1;
  nrat[0] = nrat[2] = 1e-6;
  nrat[1] = nrat[3] = 0;

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    i = (int) (rnd0() * n);
    rvn_rnddisp(xi, x[i], amp[sys]);
    dg_copy(ng, g);
    for (j = 0; j < n; j++) {
      if (j == i) continue;
      if (rvn_dist2(xi, x[j]) >= 1)
        dg_unlink(ng, i, j);
      else
        dg_link(ng, i, j);
    }

    if ( lookup ) { /* cheap lookup verseion */
      unqid_t gmapid = dg_getmapid(ng);
      nbc = dg_biconnected_lookuplow(n, gmapid);
      /* the weight `fb' is computed from the lookup table
       * for a small n, so it can be updated in every step */
      nfb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
      nnz = (nfb != 0);
    } else { /* direct and expensive computation, large n */
      nbc = dg_biconnected(ng);
      /* since computing fb is expensive,
       * we only try to find if a clique separator,
       * if so, fb == 0, otherwise fb is likely nonzero
       * so this is a good approximation of the ensemble */
      nnz = nbc ? (dg_cliquesep(ng) == 0) : 0;
    }
   
    acc1 = sys ? nnz : nbc;
    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      nz = nnz;
      if ( lookup ) fb = nfb;
    }

    /* change `sys' */
    if (sys == 0) { /* switch `sys' from 0 to 1 */
      if ( nz ) sys = 1;
    } else { /* switch `sys' from 1 to 0 */
      if ( rnd0() < rat10 ) sys = 0;
    }

    if (eql && t >= nequil) {
      printf("equilibrated at t %g, nstfb %d, lookup %d\n", t, nstfb, lookup);
      t = 0;
      it = 0;
      eql = 0;
      continue;
    }
    
    /* accumulate data after equilibration */
    if (it % nstfb == 0) {
      if ( !lookup ) fb = dg_hsfbrjw(g);
      cacc[sys] += acc1;
      hist[sys] += 1;
      fbsum[sys] += fb;
      fbnz[sys] += nz;
      wtot[sys] += 1;
    }

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      /* remove center of mass */
      rvn_rmcom(x, n);

      if (sys == 0) {
        if (rnd0() < 0.5) {
          /* the rate of biconnectivity after adding a random vertex */
          nrat[0] += 1;
          nrat[1] += rvn_voladd(x, n, xi);
        } else {
          /* the rate of biconnectivity after removing a random vertex */
          nrat[2] += 1;
          nrat[3] += dgmc_nremove(g, sg, n, &i);
        }
      }
    }

#define REPORT() \
      rathis = hist[1]/hist[0]; \
      fbav = (fbsum[0] + fbsum[1]) / (1 + fbnz[1]/fbnz[0]) / hist[0]; \
      fprintf(stderr, "D %d, n %d,%11g steps, fb%+8d(%+.8f/%+.8f|%+.8f),\n" \
          "    cacc %.4f/%.4f, nz %.4f, his1/0 %.4f, nmv+: %.8f, -: %.8f\n", \
        D, n, t, \
        fb, fbsum[0]/wtot[0], fbsum[1]/wtot[1] * rathis * rat10, fbav, \
        cacc[0]/hist[0], cacc[1]/hist[1], \
        fbnz[0]/hist[0], rathis, nrat[1]/nrat[0], nrat[3]/nrat[2]);

    if (it % nstrep == 0) {
      it = 0;
      REPORT();
    }
  }
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);

  REPORT();
  return fbav;
}



int main(int argc, char **argv)
{
  double fbav[2], nrat1[4], nrat2[4], nr1, nr2, rvir;

  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g/%g, amp %g/%g, nstfb %d, code %d-bit\n",
      D, n, (double) nsteps, nsteps2, mcamp[0], mcamp[1],
      nstfb, (int) sizeof(code_t) * 8);

  fbav[0] = mcrun(n, nequil, nsteps, mcamp,
      nstfb, lookup, nstnmv, nrat1);
  fbav[1] = mcrun(n - 1, nequil, nsteps2, mcamp,
      nstfb, lookup, nstnmv, nrat2);
  nr1 = nrat1[3]/nrat1[2];
  nr2 = nrat2[1]/nrat2[0];
  rvir = fbav[0]/fbav[1] * nr2/nr1 * (1. - n)/(2. - n)/n;
  printf("D %d, n %d, fbav %g/%g, nrat %g/%g = %g, vir%d/vir%d %g, %g, %g\n",
      D, n, fbav[0], fbav[1], nr2, nr1, nr2/nr1,
      n, n - 1, rvir, rvir * 2, rvir * ndvol(D));

  mtsave(NULL);
  return 0;
}

