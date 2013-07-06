#include <time.h>
#include "dg.h"



/* check if diagram is biconnected 
 * faster implementation */
INLINE int dg_biconnected_lookup(const dg_t *g)
{
  int k, n = g->n;
  code_t c;
  dgmap_t *m = dgmap_ + n;
  /* biconnectivity of unique diagrams */
  static int *bc[DGMAP_MAX + 1] = {NULL};

  if (bc[n] == NULL) {
    dg_t *g1 = dg_open(n);

    dgmap_init(m, n); /* compute unique maps */
    xnew(bc[n], m->ng);
    for (k = 0; k < m->ng; k++) {
      c = m->first[k];
      dg_decode(g1, &c);
      bc[n][k] = dg_biconnected(g1);
    }
    dg_close(g1);
  }
  dg_encode(g, &c);
  k = m->map[c];
  return bc[n][k];
}



static void foo(int n, int nsteps)
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
  printf("biconnectivity, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  
  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr * rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    sum += dg_biconnected_lookup(g);
  }
  printf("biconnectivity_lookup, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



int main(void)
{
  dg_t *g;
  code_t c[2];

  g = dg_open(10);
  dg_full(g);
  dg_encode(g, c);
  dg_print(g);
  printf("connected %d, biconnected %d, %#x %#x\n",
      dg_connected(g), dg_biconnected(g), c[0], c[1]);
  dg_close(g);
  foo(6, 10000000);
  return 0;
}

