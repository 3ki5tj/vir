/* iterate over permutations of n */
#include <stdio.h>
#include "dgring.h"



typedef struct {
  int npr;
  int id[64][2];
  int nr; /* the correct nr */
} edges_t;



/* test the ring contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j, err;
  double nr, nr1;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_empty(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_link(g, ref[i].id[j][0], ref[i].id[j][1]);
    if (!dg_biconnected(g))
      continue;
    nr = dg_nring_direct(g);
    nr1 = dg_nring_spec0(g, NULL, NULL, &err);
    if (fabs(nr - ref[i].nr) > 1e-3
      || (err == 0 && fabs(nr - nr1) > 1e-3)) {
      printf("n %d: model %d nr mismatch %g,%g vs %d (ref)\n",
          n, i, nr, nr1, ref[i].nr);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, # of subrings of %d reference diagrams verified\n", n, i);
}



int main(int argc, char **argv)
{
  edges_t ref4[] = {
    {5, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, 1},
    {4, {{0, 1}, {1, 3}, {2, 3}, {2, 0}}, 1},
    {4, {{0, 2}, {2, 1}, {1, 3}, {3, 0}}, 1},
    {5, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}}, 1},
    {6, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 3},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref5[] = {
    {10, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 12},
    {-1, {{0, 0}}, 0},
  };
  cmpref(4, ref4);
  cmpref(5, ref5);

  return 0;
}
