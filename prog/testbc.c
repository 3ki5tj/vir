#include <time.h>
#include "dgmap.h"



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
    bc1 = dg_biconnected_std(g);
    if (cn != ref[i].cn || bc != ref[i].bc || bc1 != ref[i].bc) {
      printf("n %d: model %d, connected %d vs %d, biconnected %d, %d vs %d (ref)\n",
          n, i, cn, ref[i].cn, bc, bc1, ref[i].bc);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, hard-sphere weights of %d reference diagrams verified\n", n, i);
}



#define LISTSIZE 100000

static void speed_connected(int n, int nsamp)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i = 0, j = 1, isamp = 0, sum = 0;
  double tsum = 0;
  dg_t *gls[LISTSIZE] = {NULL};
  int lscnt = 0, kk;

  for (t = 0; t < LISTSIZE; t++)
    gls[t] = dg_open(n);
  dg_full(g);
  for (t = 0; ; t++) {
    parsepairindex((int) (npr *rnd0()), n, &i, &j);
    if (dg_linked(g, i, j)) dg_unlink(g, i, j);
    else dg_link(g, i, j);
    dg_copy(gls[lscnt], g);
    lscnt++;
    if (t % LISTSIZE != 0) continue;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum += dg_connected(gls[kk]);
    tsum += clock() - t0;
    isamp += lscnt;
    lscnt = 0;
    if (isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("connected, n %d: time used: %gs/%d = %gmcs\n",
      n, tsum, nsamp, tsum/nsamp*1e6);
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



static void speed_biconnected(int n, int nsamp, int nedmax)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, npr = n*(n-1)/2, i = 0, j = 1, isamp = 0, ned = 0;
  int sum[3] = {0};
  double tsum[3] = {0};
  dg_t *gls[LISTSIZE] = {NULL};
  int lscnt = 0, kk;

  for (t = 0; t < LISTSIZE; t++)
    gls[t] = dg_open(n);
  dg_empty(g);
  for (i = 0; i < n; i++) dg_link(g, i, (i+1)%n);
  ned = dg_nedges(g);
  for (t = 0; ; t++) {
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
    dg_copy(gls[lscnt], g);
    lscnt++;
    if (t % LISTSIZE != 0) continue;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[0] += dg_biconnected_simple(gls[kk]);
    tsum[0] += clock() - t0;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[1] += dg_biconnected_std(gls[kk]);
    tsum[1] += clock() - t0;

#ifdef DGMAP_EXISTS
    if (n <= DGMAP_NMAX) {
      t0 = clock();
      for (kk = 0; kk < lscnt; kk++)
        sum[2] += dg_biconnected_lookup(gls[kk]);
      tsum[2] += clock() - t0;
    }
#endif

    isamp += lscnt;
    lscnt = 0;
#ifdef TESTRING
    if (ned > n) {
      dg_unlink(g, i, j);
      ned--;
    }
#endif
    if (isamp >= nsamp) break;
  }
  tsum[0] /= 1. * CLOCKS_PER_SEC * nsamp;
  tsum[1] /= 1. * CLOCKS_PER_SEC * nsamp;
  tsum[2] /= 1. * CLOCKS_PER_SEC * nsamp;
  printf("biconnected, n %d: time of %d samples simple/std/lookup: "
      "%g/%g/%gmcs, av %g/%g/%g\n",
      n, nsamp, tsum[0]*1e6, tsum[1]*1e6, tsum[2]*1e6,
      1.*sum[0]/nsamp, 1.*sum[1]/nsamp, 1.*sum[2]/nsamp);
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
  int n = DG_NMAX, nsteps = 100000;
  int n1 = 9, nsteps1 = 100000000, method = 0, nedmax = 100000;

#ifndef N
  cmpref(3, ref3);
  cmpref(4, ref4);
  cmpref(5, ref5);
#else
  n = n1 = N;
#endif

  if (argc > 1) n1 = n = atoi(argv[1]);
  if (argc > 2) nsteps1 = nsteps = atoi(argv[2]);
  printf("verification: n %d, nsteps %d\n", n, nsteps);
  verify_biconnected(n, nsteps);

  if (argc > 3) n1 = atoi(argv[3]);
  if (argc > 4) nsteps1 = atoi(argv[4]);
  if (argc > 5) nedmax = atoi(argv[5]);
  printf("speed test: n %d, nsteps %d, nedmax %d\n",
      n1, nsteps1, nedmax);
  /* n = 9, 0.0mcs on T60 */
  speed_connected(n1, nsteps1);
  /* time on T60 about 0.26mcs for n = 9 */
  speed_biconnected(n1, nsteps1, nedmax);
  return 0;
}
