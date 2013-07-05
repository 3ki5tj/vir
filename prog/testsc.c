#include <time.h>
#include "dgsc.h"


typedef struct {
  int npr;
  int id[64][2];
  int sc; /* expected sc */
} edges_t;


/* test the star contents of diagrams against the reference values */
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



static void testspeed(int n, int nsteps)
{
  int t, ipr, npr = n * (n - 1)/2, i, j, sum = 0;
  dg_t *g;
  clock_t t0;

  t0 = clock();
  g = dg_open(n);
  dg_full(g);
  for (t = 0; t < nsteps; t++) {
    /* randomly switch an edge */
    ipr = (int) (npr * rnd0());
    parsepairindex(ipr, n, &i, &j);
    if (dg_linked(g, i, j)) {
      dg_unlink(g, i, j);
      if (!dg_biconnected(g))
        dg_link(g, i, j);
    } else {
      dg_link(g, i, j);
    }
    sum += dg_rhsc(g);
  }
  printf("star content, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
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
  testspeed(5, 10000000);
  testspeed(6, 10000);
  return 0;
}
