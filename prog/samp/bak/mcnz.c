/* sampling with two ensembles
 * one ensemble samples all biconnected diagrams
 * the other ensemble samples all biconnected diagrams with nonzero `fb'
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcnz.c -lm */



#define ZCOM_PICK
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"


int n = 7; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000;
real mcamp[2] = {1.5f, 0.6f};
int nstrep = 100000000; /* interval of reporting */
double rat10 = 0.03; /* estimated rate of nonzero biconnected diagrams */
double ratcr = 0.1; /* rate of coordinate replacement */


/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-a", "%r", &mcamp[0], "MC amplitude for biconnected diagrams");
  argopt_add(ao, "-A", "%r", &mcamp[1], "MC amplitude for nonzero biconnected diagrams");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-Z", "%lf", &rat10, "estimated ratio of nonzero biconnected diagrams");
  argopt_add(ao, "-R", "%lf", &ratcr, "rate of coordinate replacement");
  argopt_parse(ao, argc, argv);
  argopt_dump(ao);
  argopt_close(ao);
  mcamp[0] /= D;
  mcamp[1] /= D;
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps, double amp[2])
{
  rvn_t *x, *nx, xi;
  int i, j, it, fb = 0, nfb = 0, acc1, sys = 0, eql;
  double fbav, fbsum[2] = {0}, nzsm = 0, t;
  dg_t *g, *ng;
  double hist[2] = {1e-6, 1e-6};
  double cacc[2] = {0}, ctot[2] = {0};
  double racc = 0, rtot = 1e-6;

  xnew(x, n + 1);
  xnew(nx, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.5 / sqrt(D)), (real) (0.5 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb(g);
  sys = 0;
  eql = 1;

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    if (sys == 0 && rnd0() < ratcr) { /* do replacement */
      if ( grepl(x, nx, g, ng) ) {
        fb = dg_hsfb(g);
        racc += 1;
      }
      rtot += 1;
    } else {
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
      acc1 = 0;
      if (sys) {
        acc1 = (nfb != 0);
      } else {
        acc1 = dg_biconnected( ng );
      }

      if ( acc1 ) { /* accept the move */
        rvn_copy(x[i], xi);
        dg_copy(g, ng);
        fb = nfb;
      }
      cacc[sys] += acc1;
      ctot[sys] += 1;
    }

    /* change `sys' */
    if (sys == 0) { /* try to turn sys from 0 to 1 */
      if (fb != 0) sys = 1;
    } else { /* try to turn sys from 1 to 0 */
      die_if (fb == 0, "fb %d must be nonzero if sys %d\n", fb, sys);
      if (rnd0() < rat10) sys = 0;
    }

    if (eql) {
      if (t >= nequil) {
        printf("equilibrated at t %g\n", t);
        t = 0;
        it = 0;
        eql = 0;
      }
      continue;
    }

    /* accumulate data after equilibration */
    hist[sys] += 1;
    fbsum[sys] += fb;
    if (sys == 0) nzsm += (fb != 0);
  }
  free(x);
  free(nx);
  dg_close(g);
  dg_close(ng);

  fbav = (fbsum[0] + fbsum[1]) / (1 + hist[1]/nzsm) / hist[0];
  fprintf(stderr, "D %d, n %d, cacc %.6f/%.6f, fb %.8f/%.8f|%.8f, "
      "nz %.6f, hisrat 1/0: %g, racc %.6f\n",
      D, n, cacc[0]/ctot[0], cacc[1]/ctot[1],
      fbsum[0]/hist[0], fbsum[1]/hist[0] * rat10, fbav,
      nzsm/hist[0], hist[1]/hist[0], racc/rtot);
  return 0;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g, amp %g/%g, code %d-bit\n",
      D, n, (double) nsteps, mcamp[0], mcamp[1], (int) sizeof(code_t) * 8);
  mcrun(n, nequil, nsteps, mcamp);
  mtsave(NULL);
  return 0;
}

