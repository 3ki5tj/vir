#include <time.h>
#include "dg.h"
#include "dgcsep.h"


/* test the RTL example */
static void testrtl(void)
{
  dg_t *g = dg_open(9), *f = dg_open(9);
  int a[9], p[9], i;

  /* RTL paper index = 9 - i */
  dg_link(g, 0, 1);
  dg_link(g, 0, 2);
  dg_link(g, 0, 3);
  dg_link(g, 1, 5);
  dg_link(g, 1, 4);
  dg_link(g, 2, 3);
  dg_link(g, 2, 6);
  dg_link(g, 3, 5);
  dg_link(g, 3, 7);
  dg_link(g, 4, 5);
  dg_link(g, 4, 8);
  dg_link(g, 5, 7);
  dg_link(g, 6, 7);
  dg_link(g, 7, 8);
  dg_print(g);
  dg_minimalorder(g, f, a, p);
  dg_print(f);
  for (i = 0; i < 9; i++) {
    printf("%d %d %d\n", i, a[i], p[i]);
  }
  dg_close(g);
  dg_close(f);
}



static void testcsep(void)
{
  dg_t *g;
  code_t c, bw;

  g = dg_open(4);
  dg_full(g);
  dg_unlink(g, 0, 2);
  dg_print(g);
  printf("clique separator: ");
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n");
  dg_close(g);

  g = dg_open(5);
  dg_full(g);
  dg_unlink(g, 0, 1);
  dg_unlink(g, 2, 3);
  dg_unlink(g, 3, 4);
  dg_print(g);
  printf("clique separator: ");
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n");
  dg_close(g);
}



static void speed_minimalorder(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i, j;

  dg_full(g);
  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    dg_minimalorder(g, NULL, NULL, NULL);
  }
  printf("dg_minimalorder, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



int main(void)
{
  testrtl();
  testcsep();
  speed_minimalorder(8, 1000000);
  return 0;
}

