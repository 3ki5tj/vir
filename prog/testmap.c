#include "dgmap.h"
#include "testutil.h"



#ifdef DGMAP_EXISTS
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
#endif



int main(void)
{
  dg_t *g;
  dgword_t c[2];

#ifdef __INTEL_COMPILER
  {
    printf("first bits: 0x14 %d, 0x80000000 %d, 0 %d\n",
        _bit_scan_forward(0x14),
        _bit_scan_forward(0x80000000),
        _bit_scan_forward(0));
  }
#endif

#ifdef DGMAP_EXISTS
  testperm(3);
#endif

  g = dg_open(10);
  dg_full(g);
  dg_encode(g, c);
  dg_print(g);
  printf("connected %d, biconnected %d, %#x %#x\n",
      dg_connected(g), dg_biconnected(g), c[0], c[1]);
  dg_close(g);

  DG_FREEMEMORIES()
  return 0;
}

