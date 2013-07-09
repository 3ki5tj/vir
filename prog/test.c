#include <time.h>
#include "dg.h"



/* test the permutation routine */
static void testperm(int n)
{
  int *p, np, i, j;

  np = dgmap_getperm(n, &p);
  for (i = 0; i < np; i++) {
    printf("%4d: ", i);
    for (j = 0; j < n; j++) printf("%d ", p[n*i+j]);
    printf("\n");
  }
  free(p);
}



static void speed_biconnected(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i, j, sum = 0;

  dg_full(g);
  dg_biconnected_lookup(g);

  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    sum += dg_biconnected(g);
  }
  printf("biconnected, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  
  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr * rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    sum += dg_biconnected_lookup(g);
  }
  printf("biconnected_lookup, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



int main(void)
{
  dg_t *g;
  code_t c[2];

  testperm(3);

  g = dg_open(10);
  dg_full(g);
  dg_encode(g, c);
  dg_print(g);
  printf("connected %d, biconnected %d, %#x %#x\n",
      dg_connected(g), dg_biconnected(g), c[0], c[1]);
  dg_close(g);
  speed_biconnected(7, 10000000);
  return 0;
}

