/* #define D to define the dimension */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"


typedef unsigned long long ulong_t;
typedef long long long_t;

int n = 4; /* order */
long_t nstepsmc = 100000000ll;
long_t nstepsacc = 0ll;



/* compute the sign of the virial coefficient */
static double mcrun(int n, long_t nsteps, double *ring)
{
  rvn_t *x, xi;
  int par, i, j, ipr, npr = n * (n - 1)/2, nedg;
  long_t cnt[2] = {0, 0}, nring = 0, tot, t;
  long_t gtot = 0, grej = 0;
  long_t ctot = 0, cacc = 0;
  real amp = (real) 2.0/D, rate;
  dg_t *g;

  xnew(x, n);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-.5 / sqrt(D)), (real) (.5 / sqrt(D)) );
  g = dg_open(n);
  dg_full(g); /* fully connected diagram */
  nedg = npr;
  par = npr % 2;

  for (t = 0; t < nsteps; t++) {
    if (rnd0() < 0.1) /* change the underlying coordinates */
    {
      i = 1 + (int) (rnd0() * (n - 1));
      rvn_rnddisp(xi, x[i], amp);
      /* check if all constraints are satisfied */
      for (j = 0; j < n; j++)
        if (dg_linked(g, i, j) && rvn_dist2(xi, x[j]) > 1)
          break;
      ctot++;
      if (j == n) { /* compatible diagram */
        cacc++;
        rvn_copy(x[i], xi);
      }
    }
    else  /* change the graph topology */
    {
      const double prm = 1.0;

      ipr = npr * rnd0();
      parsepairindex(ipr, n, &i, &j);
      gtot++;
      if ( dg_linked(g, i, j) ) { /* try to unlink an existing graph */
        if ( prm >= 1 || rnd0() < prm ) {
          dg_unlink(g, i, j);
          if ( !dg_biconnected(g) ) {
            grej++;
            dg_link(g, i, j); /* link back */
          } else {
            par = !par;
            nedg--;
          }
        }
      } else { /* add a link */
        if ( (prm <= 1 || rnd0() < 1/prm)
          && rvn_dist2(x[i], x[j]) < 1) {
          dg_link(g, i, j);
          par = !par;
          nedg++;
        }
      }
    }
    cnt[par]++;
    if (nedg == n) nring++;
  }
  tot = cnt[0] + cnt[1];
  rate = 1.*(cnt[0] - cnt[1])/tot;
  *ring = 1.*nring/tot;
  printf("even %lld, odd %lld, del %g, gacc %g, cacc %g, nring %g\n",
      cnt[0], cnt[1], rate,
      1. - 1.*grej/gtot, 1.*cacc/ctot, *ring);
  free(x);
  dg_close(g);
  return rate;
}


/* compute the acceptence rate of the ring diagram */
static double calcringacc(int n, long_t nsteps)
{
  long_t t, acc = 0;
  int i;
  rvn_t x, u;

  for (t = 0; t < nsteps; t++) {
    rvn_zero(x);
    for (i = 0; i < n - 1; i++)
      rvn_inc(x, rvn_rndball0(u));
    if (rvn_sqr(x) < 1) acc++;
  }
  return 1. * acc / nsteps;
}



int main(int argc, char **argv)
{
  double rate, ringr, ringacc, vir;

  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nstepsmc = (long_t) atof(argv[2]);
  if (argc > 3) nstepsacc = (long_t) atof(argv[3]);
  if (nstepsacc <= 0) nstepsacc = nstepsmc / 100;
  printf("n = %d, D = %d, nsteps = %lld, nstepsacc = %lld\n",
      n, D, nstepsmc, nstepsacc);

  rate = mcrun(n, nstepsmc, &ringr);
  ringacc = calcringacc(n, nstepsacc);

  vir = rate * ringacc / ringr;
  vir *= (1. - n) / (2. * n);
  printf("rate %.6f, ringr %.6f, ringacc %.6f, vir %.6f, %.6f, %g\n",
      rate, ringr, ringacc, vir,
      vir * pow(2, n - 1), vir * pow(8, n - 1));
  mtsave(NULL);
  return 0;
}

