#include "dgcryr.h"


static void test1(void)
{
  dg_t *g;
  int i, j, n = 4;

  g = dg_open(n);
  dg_link(g, 0, 1);
  dg_link(g, 0, 2);
  dg_link(g, 0, 3);
  dg_link(g, 1, 2);
  dg_link(g, 1, 3);
  dg_print(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ )
      printf("i %d, j %d, %g\n", i, j, dgsc_yiter(g, NULL, i, j));
  dg_close(g);
}



int main(void)
{
  test1();
  return 0;
}
