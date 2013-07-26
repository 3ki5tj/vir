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
  return 0;
}

