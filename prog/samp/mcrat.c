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
#include "dgmcutil.h"



int n = 4; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000;
real mcamp = 2.f;
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 20; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp, "MC amplitude");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-z", "%d", &nstnmv, "interval of the n-move");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
  if (nstfb <= 0) /* frequency of computing the weight */
    nstfb = (n > DGMAP_NMAX) ? 100 : 1;
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps, real amp,
    int nstfb, int nstnmv, int plus, double *nrat)
{
  rvn_t *x, xi;
  int i, j, it, fb = 0;
  double fbav, nzsum = 0, fbsum = 0, t, acc = 0, tot = 0, nsum = 0, ntot = 0;
  dg_t *g, *ng, *sg;

  xnew(x, n + 1);
  amp *= (real) 1.0/D;
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb(g);

  for (it = 1, t = 1; t <= nsteps + nequil; t += 1, it++) {
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
      /* the weight `fb' is computed from the lookup table
       * for a small n, so it can be updated in every step */
      if (n <= DGMAP_NMAX)
        fb = dg_hsfb(g);
      acc += 1.;
    }
   
    if (t < nequil) {
      continue;
    } else if (t == nequil) {
      printf("equilibrated at t %g\n", t);
    }
    
    /* accumulate data after equilibration */
    if (n <= DGMAP_NMAX) {
#define ACCUM() { \
      tot += 1.;  \
      fbsum += fb;  \
      nzsum += (fb != 0); }
      ACCUM();
    } else if (it % nstfb == 0) {
      fb = dg_hsfb(g);
      ACCUM();
    }

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      /* remove center of mass */
      rvn_rmcom(x, n);

      if (plus) {
        ntot += 1;
        nsum += rvn_voladd(x, n, xi);
      } else {
        /* compute the probability of biconnectivity after removing a random vertex */
        ntot += 1;
        nsum += dgmc_nremove(g, sg, n, &i);
      }

      if (it % nstrep == 0) {
        it = 0;
        fprintf(stderr, "n %d, %g steps, tot %g, acc %.8f, "
            "fb %d(%.8f), nz %.8f, nc %.8f\n",
          n, t, tot, acc/t, fb, fbsum/tot, nzsum/tot, nsum/ntot);
      }
    }
  }
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);

  fbav = fbsum / tot;
  *nrat = nsum / ntot;
  fprintf(stderr, "n %d, acc %.8f, fb %.8f, nz %g, nc %.8f\n",
      n, acc/nsteps, fbav, nzsum/tot, *nrat);
  return fbav;
}



int main(int argc, char **argv)
{
  double fbav[2], nrat[2], rvir;

  doargs(argc, argv);
  printf("n %d, D %d, nsteps %g, amp %g, nstfb %d, code %d-bit\n",
      n, D, (double) nsteps, mcamp, nstfb, (int) sizeof(code_t) * 8);

  fbav[0] = mcrun(n, nequil, nsteps, mcamp,
      nstfb, nstnmv, 0, nrat);
  /* the smaller system converges faster, so we spend less time on it */
  fbav[1] = mcrun(n - 1, nequil, nsteps / 4, mcamp,
      nstfb, nstnmv, 1, nrat + 1);
  rvir = fbav[0]/fbav[1] * (nrat[1]/nrat[0]) * (1. - n)/(2. - n)/n;
  printf("n %d, fbav %g/%g, nrat %g/%g = %g, vir%d/vir%d %g, %g, %g\n",
      n, fbav[0], fbav[1], nrat[0], nrat[1], nrat[0]/nrat[1],
      n, n - 1, rvir, rvir * 2, rvir * ndvol(D));

  mtsave(NULL);
  return 0;
}

