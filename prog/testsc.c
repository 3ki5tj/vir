#include "dgmap.h"
#include "dgsc.h"
#include "dgrjw.h"
#include "testutil.h"



typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int sc; /* the correct sc */
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j;
  double sc;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    if (!dg_biconnected(g))
      continue;
    sc = dgsc_do(g);
    if (fabs(sc - ref[i].sc) > 0.001) {
      printf("n %d: model %d sc mismatch %g vs %d (ref)\n",
          n, i, sc, ref[i].sc);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, star contents of %d reference diagrams verified\n", n, i);
}



/* test the speed of computing star content
 * see the function with the same name for a more advanced version */
static void testspeed(int n, int nsamp, int nedmax, char method)
{
  int t, ned, nequil = 1000, eql = 1, isamp = 0, good = 0, tot = 0;
  dg_t *g;
  clock_t t0;
  double sum, tsum = 0, rnp = 0.1;

  g = dg_open(n);
  dg_full(g);
#ifdef DGMAP_EXISTS
  if (method == 'l' && n < DGMAP_NMAX) {
    int ig;
    double sc1, sc2;
    dg_t *g2;

    t0 = clock();
    dg_sc(g); /* automatically activate the look up table */
    printf("star content, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
    /* compare with the RJW result */
    g2 = dg_open(g->n);
    for (ig = 0; ig < dgmap_[n].ng; ig++) {
      dgword_t code = dgmap_[n].first[ig];
      dg_decode(g2, &code);
      sc1 = dg_sc(g2);
      sc2 = dg_fb(g2);
      if (dg_nedges(g2) % 2 == 1) sc2 *= -1;
      if (fabs(sc1 - sc2) > 1e-3) {
        printf("sc1 %g, sc2 %g\n", sc1, sc2);
        dg_print(g2);
        exit(1);
      }
    }
    dg_close(g2);
  }
#endif

  ned = dg_nedges(g);
  for (t = 0; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge(g, &ned, nedmax, rnp);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */

    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if ( ned > nedmax ) continue;
    if ( dg_cliquesep(g) ) continue;

    t0 = clock();
    sum += dg_sc(g);
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("star content, n %d, method %c, sum %g, time used: %gs/%d = %gms\n",
      n, method, sum, tsum, nsamp, tsum/nsamp*1000);
  dg_close(g);
}



int main(void)
{
  edges_t ref5[] = { /* edges are wiggly lines */
    {0, {{0, 0}}, -6},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, 3},
    {2, {{0, 1}, {1, 3}}, 0},
    {3, {{0, 1}, {2, 3}, {3, 4}}, -2},
    {3, {{0, 1}, {1, 3}, {2, 4}}, -2},
    {3, {{0, 1}, {1, 3}, {1, 4}}, 0},
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
  int nedmax = 1000;

  cmpref(5, ref5);
  cmpref(6, ref6);
  testspeed(5, 100000,  nedmax, 'l');
  testspeed(6, 100000,  nedmax, 'l');
  testspeed(7, 100000,  nedmax, 'l');
  testspeed(8, 1000000, nedmax, 'l');
  testspeed(8, 1000,    nedmax, 'd');
  return 0;
}
