#include "dg.h"

int main(void)
{
  dg_t *g;

  g = dg_open(5);
  dg_link(g, 0, 2);
  dg_link(g, 0, 3);
  dg_link(g, 0, 4);
  dg_link(g, 1, 2);
  dg_link(g, 1, 3);
  dg_link(g, 1, 4);
  dg_print(g);
  printf("connected %d, biconnected %d\n", dg_connected(g), dg_biconnected(g));
  dg_close(g);
  return 0;
}

