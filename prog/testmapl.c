#include "dgmapl.h"


static void foo(int k, int n, int nsteps)
{
  dg_t *g, *ng;
  int i, j, t, err = 0, st[DGMAPL_NMAX + 1];
  clock_t t0;
  double tsum = 0;

  g = dg_open(n);
  ng = dg_open(n);
  dg_full(g);
  for (t = 0; t < nsteps; t++) {
    i = randpair(n, &j);
    if ( dg_linked(g, i, j) ) {
      dg_unlink(g, i, j);
      if (!dg_biconnected(g))
        dg_link(g, i, j); /* link back */
    } else {
      dg_link(g, i, j);
    }
    t0 = clock();
    err += (dgmapl_getchain(g, k, st) == 0);
    tsum += clock() - t0;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("getchain() err %d,%g, %g/%d = %gmcs\n",
    err, 1.*err/nsteps, tsum, nsteps, tsum/nsteps*1e6);
  dg_close(g);
  dg_close(ng);
}


int main(int argc, char **argv)
{
  int n = 9, nsteps = 1000000, k = 4;
  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsteps = atoi(argv[2]);
  if (argc >= 4) k = atoi(argv[3]);
  foo(k, n, nsteps);
  return 0;
}

