#include "dg.h"

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
  return 0;
}

