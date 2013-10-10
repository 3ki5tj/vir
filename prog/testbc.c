#include <time.h>
#include "dg.h"


typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int cn, bc;
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j, cn, bc, bc1;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_empty(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_link(g, ref[i].id[j][0], ref[i].id[j][1]);
    cn = dg_connected(g);
    bc = dg_biconnected(g);
    bc1 = dg_biconnected(g);
    if (cn != ref[i].cn || bc != ref[i].bc || bc != ref[i].bc) {
      printf("n %d: model %d, connected %d vs %d, biconnected %d, %d vs %d (ref)\n",
          n, i, cn, ref[i].cn, bc, bc1, ref[i].bc);
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
  int t, npr = n*(n-1)/2, i = 0, j = 1, sum = 0;
  double tsum = 0;

  dg_full(g);
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    t0 = clock();
    sum += dg_connected(g);
    tsum += clock() - t0;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("connected, n %d: time used: %gs/%d = %gmcs\n",
      n, tsum, nsteps, tsum/nsteps*1e6);
  dg_close(g);
}



/* check biconnectivity by two algorithms */
static void verify_biconnected(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  int t, npr = n*(n-1)/2, i = 0, j = 1, bc1, bc2;

  dg_full(g);
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    bc1 = dg_biconnected(g);
    bc2 = dg_biconnected_std(g);
    if (bc1 != bc2) {
      printf("bc1 %d vs bc2 %d\n", bc1, bc2);
      dg_print(g);
      exit(1);
    }
  }
  printf("verified biconnectivity of %d diagrams\n", nsteps);
  dg_close(g);
}



static void speed_biconnected(int n, int nsteps,
    int method, int nedmax)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i = 0, j = 1, sum = 0, ned = 0;
  double tsum = 0;

  die_if (method == 'l' && n > DGMAP_NMAX,
    "n %d cannot use the lookup method\n", n);

  dg_empty(g);
  for (i = 0; i < n; i++) dg_link(g, i, (i+1)%n);
  ned = dg_nedges(g);
  for (t = 0; t < nsteps; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
#ifdef TESTRING /* test ring-like conformations */
    if ( !dg_linked(g, i, j) ) {
      dg_link(g, i, j);
      ned++;
    }
#else
    if (dg_linked(g, i, j)) {
      dg_unlink(g, i, j);
      ned--;
    } else if (ned < nedmax) {
      dg_link(g, i, j);
      ned++;
    }
#endif
    t0 = clock();
    if (method == 's')
      sum += dg_biconnected_std(g);
    else if (method == 'l')
      sum += dg_biconnected_lookup(g);
    else
      sum += dg_biconnected(g);
    tsum += clock() - t0;
#ifdef TESTRING
    if (ned > n) {
      dg_unlink(g, i, j);
      ned--;
    }
#endif
  }
  tsum /= CLOCKS_PER_SEC;
  printf("biconnected, n %d: time used: %gs/%d = %gmcs, av %g\n",
      n, tsum, nsteps, tsum / nsteps * 1e6, 1. * sum / nsteps);
  dg_close(g);
}




int main(int argc, char **argv)
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
    {-1, {{0, 0}}, 1, 1},
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
    {6, {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}}, 1, 1},
    {-1, {{0, 0}}, 1, 1},
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
    {-1, {{0, 0}}, 1, 1},
  };
  int n = DG_NMAX, nsteps = 1000000;
  int n1 = 9, nsteps1 = 100000000, method = 0, nedmax = 100000;

  cmpref(3, ref3);
  cmpref(4, ref4);
  cmpref(5, ref5);

  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nsteps = atoi(argv[2]);
  printf("verification: n %d, nsteps %d\n", n, nsteps);
  verify_biconnected(n, nsteps);

  if (argc > 3) n1 = atoi(argv[3]);
  if (argc > 4) nsteps1 = atoi(argv[4]);
  if (argc > 5) method = argv[5][0];
  if (argc > 6) nedmax = atoi(argv[6]);
  printf("speed test: n %d, nsteps %d, method %c, nedmax %d\n",
      n1, nsteps1, method, nedmax);
  speed_connected(n1, nsteps1);
  speed_biconnected(n1, nsteps1, method, nedmax);
  return 0;
}
