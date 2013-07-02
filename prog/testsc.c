#include "dgsc.h"


typedef struct {
  int npr;
  int id[64][2];
  int sc; /* expected sc */
} edges_t;



static void test(int n, edges_t *ref)
{
  int i, j, sc;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    sc = dg_rhsc(g);
    if (!dg_biconnected(g))
      continue;
    if (sc != ref[i].sc) {
      printf("n %d: model %d sc mismatch %d vs %d\n",
          n, i, sc, ref[i].sc);
      dg_print(g);
    }
  }
  dg_close(g);
}



int main(void)
{
  edges_t ref5[] = {
    {0, {{0, 0}}, -6},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, 3},
    {3, {{0, 1}, {2, 3}, {3, 4}}, -2},
    {4, {{0, 1}, {2, 3}, {3, 4}, {2, 4}}, 1},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, 1},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref6[] = {
    {0, {{0, 0}}, 24},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, -12},
    {3, {{0, 1}, {2, 3}, {3, 4}}, 8},
    {4, {{0, 1}, {2, 3}, {3, 4}, {2, 4}}, -4},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -4},
    {-1, {{0, 0}}, 0},
  };

  test(5, ref5);
  test(6, ref6);
  return 0;
}
