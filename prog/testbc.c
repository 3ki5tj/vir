#include "dgrjw.h"

#include <time.h>
#include "dgsc.h"


typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int cn, bc;
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j, cn, bc;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_empty(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_link(g, ref[i].id[j][0], ref[i].id[j][1]);
    cn = dg_connected(g);
    bc = dg_biconnected(g);
    if (cn != ref[i].cn || bc != ref[i].bc) {
      printf("n %d: model %d, connected %d vs %d, biconnected %d vs %d (ref)\n",
          n, i, cn, ref[i].cn, bc, ref[i].bc);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, hard-sphere weights of %d reference diagrams verified\n", n, i);
}



static void speed_connected(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i, j, sum = 0;

  dg_full(g);

  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    sum += dg_connected(g);
  }
  printf("connected, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



static void speed_biconnected(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i, j, sum = 0;

  dg_full(g);

  t0 = clock();
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    sum += dg_biconnected(g);
  }
  printf("biconnected, n %d: time used: %gs/%d\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);

  if (n <= DGMAP_NMAX) {
    dg_biconnected_lookup(g);
    t0 = clock();
    for (t = 0; t < nsteps; t++) {
      parsepairindex((int) (npr * rnd0()), n, &i, &j);
      if (dg_linked(g, i, j)) dg_unlink(g, i, j);
      else dg_link(g, i, j);
      sum += dg_biconnected_lookup(g);
    }
    printf("biconnected_lookup, n %d: time used: %gs/%d\n",
        n, 1.*(clock() - t0) / CLOCKS_PER_SEC, nsteps);
  }
  dg_close(g);
}



int main(void)
{
  edges_t ref3[] = {
    {0, {{0, 0}}, 0, 0},
    {1, {{0, 1}}, 0, 0},
    {1, {{0, 2}}, 0, 0},
    {1, {{1, 2}}, 0, 0},
    {2, {{0, 1}, {1, 2}}, 1, 0},
    {2, {{0, 1}, {0, 2}}, 1, 0},
    {2, {{0, 2}, {1, 2}}, 1, 0},
    {3, {{0, 1}, {1, 2}, {0, 2}}, 1, 1},
    {-1, {{0, 0}}, 1},
  };

  edges_t ref4[] = {
    {0, {{0, 0}}, 0, 0},
    {1, {{0, 1}}, 0, 0},
    {1, {{2, 3}}, 0, 0},
    {2, {{0, 1}, {1, 2}}, 0, 0},
    {2, {{0, 2}, {3, 2}}, 0, 0},
    {2, {{0, 2}, {1, 3}}, 0, 0},
    {3, {{0, 1}, {1, 2}, {2, 3}}, 1, 0},
    {3, {{0, 1}, {1, 2}, {0, 2}}, 0, 0},
    {3, {{0, 1}, {0, 2}, {0, 3}}, 1, 0},
    {4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}}, 1, 0},
    {4, {{0, 1}, {1, 3}, {0, 3}, {1, 2}}, 1, 0},
    {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, 1, 1},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}}, 1, 1},
    {6, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 4}}, 1, 1},
    {-1, {{0, 0}}, 1},
  };

  edges_t ref5[] = {
    {0, {{0, 0}}, 0, 0},
    {1, {{0, 1}}, 0, 0},
    {1, {{2, 3}}, 0, 0},
    {2, {{0, 1}, {1, 2}}, 0, 0},
    {2, {{0, 2}, {3, 2}}, 0, 0},
    {2, {{0, 2}, {1, 3}}, 0, 0},
    {3, {{0, 1}, {1, 2}, {2, 3}}, 0, 0},
    {3, {{0, 1}, {1, 2}, {0, 2}}, 0, 0},
    {3, {{0, 1}, {0, 2}, {0, 3}}, 0, 0},
    {4, {{0, 1}, {0, 2}, {0, 3}, {1, 2}}, 0, 0},
    {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, 0, 0},
    {4, {{0, 1}, {1, 2}, {1, 3}, {1, 4}}, 1, 0},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, 1, 1},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {3, 0}}, 1, 0},
    {6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {1, 4}}, 1, 1},
    {6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 3}}, 1, 1},
    {6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 1}, {0, 3}}, 1, 1},
    {7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 3}, {1, 4}}, 1, 1},
    {7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 2}, {0, 3}}, 1, 1},
    {8, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 2}, {0, 3}, {1, 4}}, 1, 1},
    {-1, {{0, 0}}, 1},
  };

  cmpref(3, ref3);
  cmpref(4, ref4);
  cmpref(5, ref5);
  speed_connected(9, 100000000);
  speed_biconnected(9, 100000000);
  return 0;
}
