#include "dgmap.h"
#include "dgrjwb.h"
#include "testutil.h"


static void foo()
{
  int n = 4, r1, z1;
  unsigned vs, sz = MKBIT(n);
  dg_t *g;

  g = dg_open(n);
  dg_full(g);
  //dg_unlink(g, 0, 1);
  //dg_unlink(g, 2, 3);
  printf("%g\n", (double) dgrjwb_fb(g));
  for ( vs = 1; vs < sz; vs++ ) {
    dgvs_fprint0(vs, n, NULL, NULL, NULL);
    z1 = n - 1;
    for (r1 = 0; r1 <= z1; r1++)
      printf(" %4g", (double) dgrjwb_arr_[sz*n*(n%2) + vs*n + r1]);
    printf("\n");
  }
  dg_close(g);
}

int main(void)
{
  foo();
  DG_FREEMEMORIES()
  return 0;
}
