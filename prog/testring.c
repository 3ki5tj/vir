/* iterate over permutations of n */
#include "dgring.h"
#include "testutil.h"



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



static void testspeed(int n, int nsamp, int nedmax)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  double tsum = 0, sum = 0;
  double rnp = 0.1; /* rate for moves increasing edges */

  printf("speed test for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  dg_full(g);
  ned = dg_nedges(g);
  for (t = 1; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge(g, &ned, nedmax, rnp);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if ( ned > nedmax ) continue;

    t0 = clock();
    sum += dg_nring_direct(g);
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("dg_nring: n %d; ave. %g; time used: %gs/%d = %gmcs\n",
      n, sum/nsamp, tsum, nsamp, tsum/nsamp*1e6);
  dg_close(g);
}



int main(int argc, char **argv)
{
  edges_t ref4[] = {
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, 1},
    {4, {{0, 1}, {1, 3}, {2, 3}, {2, 0}}, 1},
    {4, {{0, 2}, {2, 1}, {1, 3}, {3, 0}}, 1},
    {5, {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {1, 2}, {1, 3}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {2, 3}}, 1},
    {5, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}}, 1},
    {6, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 3},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref5[] = {
    {10, {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}}, 12},
    {-1, {{0, 0}}, 0},
  };
  int n = 9, nsamp = 100000, nedmax = 1000000;

#ifdef N
  n = N;
#endif

#ifndef N
  cmpref(4, ref4);
  cmpref(5, ref5);
#endif

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  /* with default setting, nedmax = inf, N not predefined
   * T60 timing Oct. 30 2013
   * 2.7mcs, n = 7
   * 7.3mcs, n = 8
   * 25mcs,  n = 9
   * 103mcs, n = 10
   * 487mcs, n = 11 */
  testspeed(n, nsamp, nedmax);
  return 0;
}
