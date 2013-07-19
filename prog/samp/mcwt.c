/* sampling with several ensembles
 * ensemble 0 samples all biconnected diagrams
 * ensemble 1 samples all biconnected diagrams with nonzero `fb'
 * ensemble 2 samples all nonzero  biconnected diagrams with the weight 1/|fb|
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcnz.c -lm */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"

#define NSYS 3


int n = 5; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000;
real mcamp[NSYS] = {1.5f, 0.6f, 0.6f};
int nstfb = 0; /* interval of evaluting the weight */
int nstrep = 100000000; /* interval of reporting */
double Z[NSYS] = {1, 0.3, 0.6};

/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int i;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp[0], "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-A", "%r", &mcamp[1], "MC amplitude for nonzero biconnected diagrams");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
  if (nstfb <= 0) /* frequency of computing the weight */
    nstfb = (n > DGMAP_NMAX) ? 100 : 1;
  for (i = 2; i < NSYS; i++) mcamp[i] = mcamp[1];
  for (i = 0; i < NSYS; i++) mcamp[i] /= D;
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps,
    double amp[], int nstfb)
{
  rvn_t *x, xi;
  int i, j, it, fb = 0, nfb = 0, acc1, sys = 0, eql, wt, nwt;
  /* probability of fb being nonzero, and that of non-clique center */
  double t, nzsum = 0;
  dg_t *g, *ng;
  double fbsum[NSYS] = {0}, fbnz = 0, fbav;
  double tsum1 = 0, tsum2 = 0, tacc1 = 0, tacc2 = 0;
  double acc[NSYS] = {0}, hist[NSYS] = {1e-6, 1e-6, 1e-6};

  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb(g);
  wt = abs(fb);
  sys = 0;
  eql = 1;

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    i = (int) (rnd0() * n);
    rvn_rnddisp(xi, x[i], (real) amp[sys]);
    dg_copy(ng, g);
    for (j = 0; j < n; j++) {
      if (j == i) continue;
      if (rvn_dist2(xi, x[j]) >= 1)
        dg_unlink(ng, i, j);
      else
        dg_link(ng, i, j);
    }
    
    nfb = (g->c[i] != ng->c[i]) ? dg_hsfb(ng) : fb;
    nwt = abs(nfb);
    acc1 = 0;
    if (sys == 0) {
      acc1 = dg_biconnected( ng );
    } else if (sys == 1) {
      acc1 = (nfb != 0);
    } else {
      if (nfb != 0) {
        if (nwt > wt) acc1 = 1;
        else acc1 = (rnd0() < 1.*nwt/wt);
      } else acc1 = 0;
    }

    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      fb = nfb;
      wt = nwt;
    }

    /* change `sys' */
    if (sys == 0) { /* try to turn sys from 0 to 1 */
      if (rnd0() < 0.5)
        if (fb != 0) sys = 1;
    } else if (sys == 1) { /* try to turn sys from 1 to 0 or 2 */
      if (rnd0() < 0.5) {
        if (rnd0() < Z[1]/Z[0]) sys = 0;
      } else {
        /* the weight in system 2 is fb/Z[2], that in system 1 is 1/Z[1] */
        tsum1 += 1;
        wt = abs(fb);
        if (rnd0() < wt*Z[1]/Z[2]) {
          sys = 2;
          tacc1 += 1;
        }
      }
    } else { /* turn sys from 2 to 1 */
      if (rnd0() < 0.5) {
        tsum2 += 1;
        wt = abs(fb);
        if (rnd0() < Z[2]/Z[1]/wt) {
          sys = 1;
          tacc2 += 1;
        }
      }
    }
     
    if (eql && t >= nequil) {
      printf("equilibrated at t %g\n", t);
      t = 0;
      it = 0;
      eql = 0;
      continue;
    }
    
    /* accumulate data after equilibration */
    acc[sys] += acc1;
    hist[sys] += 1;
    if (sys == 0 && fb != 0) fbnz += 1;
    wt = (sys == 2) ? abs(fb) : 1;
    fbsum[sys] += fb / wt;
    
#define REPORT() \
  fbav = (fbsum[0] + fbsum[1]) / (1 + hist[1]/fbnz) / hist[0]; \
  fprintf(stderr, "D %d, n %d, acc %.6f/%.6f/%.6f, fb %.8f/%.8f|%.8f %.8f\n" \
      "  t%11g nz %.6f, hisrat 1/0: %g, 2/1: %g, tacc1 %g, tacc2 %g\n", \
      D, n, acc[0]/hist[0], acc[1]/hist[1], acc[2]/hist[2], \
      fbsum[0]/hist[0], fbsum[1]/hist[0] * (Z[1]/Z[0]), fbav, \
      fbsum[2]/hist[2] * ((tacc1/tsum1)/(tacc2/tsum2)) * (Z[2]/Z[1]) * (fbnz/hist[0]), \
      t, fbnz/hist[0], hist[1]/hist[0], hist[2]/hist[1], tacc1/tsum1, tacc2/tsum2);

    if (it % nstrep == 0) {
      it = 0;
      REPORT();
    }
  }
  free(x);
  dg_close(g);
  dg_close(ng);

  REPORT();
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g, amp %g/%g, nstfb %d, code %d-bit\n",
      D, n, (double) nsteps, mcamp[0], mcamp[1], nstfb, (int) sizeof(code_t) * 8);
  mcrun(n, nequil, nsteps, mcamp, nstfb);
  mtsave(NULL);
  return 0;
}

