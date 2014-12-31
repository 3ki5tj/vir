/* accumulate statistics
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcstat.c -lm */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"



int n = 7; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000;
real mcamp = 2.f;
int nstfb = 0; /* interval of evaluting the weight */
int nstrep = 100000000; /* interval of reporting */
int fbcutoff = -1; /* |fb| cutoff */


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
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-c", "%d", &fbcutoff, "cutoff of fb");
  argopt_parse(ao, argc, argv);

  if (nstfb <= 0) /* frequency of computing the weight */
    nstfb = (n > DGMAP_NMAX) ? 100 : 10;
  if (fbcutoff < 0) {
    static int cut[] = {0, 0, 0, 0, 1, 2, 10, 10, 10};
    fbcutoff = (n > 8) ? 10 : cut[n];
  }
  argopt_dump(ao);
  argopt_close(ao);
}



/* compute the radius of gyration */
static double getrad(rvn_t *x, int n)
{
  int i, im = -1;
  double d2max = 1e9, d2;

  /* find the vertex closest to the center */
  for (i = 0; i < n; i++)
    if ((d2 = rvn_sqr(x[i])) < d2max) {
      im = i;
      d2max = d2;
    }

  for (d2 = 0, i = 0; i < n; i++) {
    if (i == im) continue;
    d2 += rvn_dist2(x[i], x[im]);
  }
  return d2 / (n - 1);
}



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps, real amp,
    int nstfb)
{
  rvn_t *x, xi;
  int i, j, it, fb = 0, nclmax = 0, ned, ncl;
  double fbav, fbsum = 0, Fbsum = 0, t, acc = 0, tot = 0;
  /* probability of fb being nonzero, and that of non-clique center */
  double nzsum = 0, ncsum = 0, fbsumx = 0, cutcnt = 0;
  double Hcl[DG_NMAX + 1] = {0}, hsm;
  double Hed[DG_NMAX * DG_NMAX + 1][3] = {{0}};
  double d2sm = 0;
  dg_t *g, *ng;

  xnew(x, n + 1);
  amp *= (real) 1.0/D;
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  mkgraph(g, x, n);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb(g);

  for (it = 1, t = 1; t <= nsteps + nequil; t += 1, it++) {
    i = (int) (rand01() * n);
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
    if (it % nstfb != 0) continue;

    tot += 1.;
    ned = dg_nedges(g);
    Hed[ned][0] += 1;
    ncl = dg_ncsep(g);
    fb = dg_hsfb(g);
    d2sm += getrad(x, n);
    if (ncl > nclmax) nclmax = ncl;
    Hcl[ncl] += 1.;
    if (ncl == 0) {
      Hed[ned][1] += 1;
      Hed[ned][2] += fb;
    }
    ncsum += (ncl == 0);
    nzsum += (fb != 0);
    fbsum += fb;
    fbsumx += (abs(fb) > fbcutoff) ? fb : 0;
    cutcnt += (abs(fb) > fbcutoff);
    Fbsum += fabs(fb);
  }
  free(x);
  dg_close(g);
  dg_close(ng);

  fbav = fbsum / tot;
  printf("n %d, acc %.8f, fb %.8f, |fb| %.8f, %.8f,\n"
      "fb with a magnitude cutoff %.8f (cutoff %d; frac %g, %g)\n"
      "nz %.8f, nc %.8f nz/nc %.8f, rad %.6f\n",
      n, acc/nsteps, fbsum/tot, Fbsum/tot, Fbsum/nzsum,
      fbsumx/tot, fbcutoff, cutcnt/tot, cutcnt/nzsum,
      nzsum/tot, ncsum/tot, nzsum/(ncsum + 1e-5), d2sm/tot);

  printf("number of clique separators:\n");
  for (hsm = 0, i = 0; i <= nclmax; i++) hsm += Hcl[i];
  for (i = 0; i <= nclmax; i++)
    printf("%2d %.6f\n", i, Hcl[i]/hsm);
  printf("\n");

  { /* print out fb as a function of the # of edges */
    double hsm2, xx, yy;

    printf("nedgs     his     hisnz          fbnz         <fb>\n");
    for (hsm = hsm2 = 0, i = n; i <= n*(n-1)/2; i++) {
      hsm += Hed[i][0]; /* biconnected diagrams */
      hsm2 += Hed[i][1]; /* non-clique-separable */
    }
    for (yy = 0, i = n; i <= n*(n-1)/2; i++) {
      xx = (Hed[i][1] > 0 ? Hed[i][2]/Hed[i][1] : 0);
      yy += xx * Hed[i][1] / hsm;
      if (Hed[i][0] > 0)
        printf("%4d  %.7f %.7f %+14.7f %+12.7f %+12.7f\n",
            i, Hed[i][0]/hsm, Hed[i][1]/hsm2, xx, xx*Hed[i][1]/hsm, yy);
    }
  }
  return fbav;
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  printf("n %d, D %d, nsteps %g, amp %g, nstfb %d, code %d-bit\n",
      n, D, (double) nsteps, mcamp, nstfb, (int) sizeof(dgword_t) * 8);
  mcrun(n, nequil, nsteps, mcamp, nstfb);
  mtsave(NULL);
  return 0;
}

