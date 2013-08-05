#include <time.h>
#include "dgsc.h"
#include "dgrjw.h"



typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int sc; /* the correct sc */
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j, sc;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    sc = dg_rhsc_low(g);
    if (!dg_biconnected(g))
      continue;
    if (sc != ref[i].sc) {
      printf("n %d: model %d sc mismatch %d vs %d (ref)\n",
          n, i, sc, ref[i].sc);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, star contents of %d reference diagrams verified\n", n, i);
}



static void testspeed(int n, int nsteps, int lookup)
{
  int t, ipr, npr = n * (n - 1)/2, i, j, sum = 0, nequil = 1000;
  dg_t *g;
  clock_t t0, t1;

  t0 = clock();
  g = dg_open(n);
  dg_full(g);
  if (lookup) {
    int ig, sc1, sc2;
    dg_t *g2;

    dg_rhsc(g); /* initialization */
    printf("star content, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
    /* compare with the RJW result */
    g2 = dg_open(g->n);
    for (ig = 0; ig < dgmap_[n].ng; ig++) {
      code_t code = dgmap_[n].first[ig];
      dg_decode(g2, &code);
      sc1 = dg_rhsc(g2);
      sc2 = dg_hsfb(g2);
      if (dg_nedges(g2) % 2 == 1) sc2 *= -1;
      if (sc1 != sc2) {
        printf("sc1 %d, sc2 %d\n", sc1, sc2);
        dg_print(g2);
        exit(1);
      }
    }
    dg_close(g2);
  }
  for (t = 0; t < nequil + nsteps; t++) {
    /* randomly switch an edge */
    ipr = (int) (npr * rnd0());
    parsepairindex(ipr, n, &i, &j);
    if (dg_linked(g, i, j)) {
      dg_unlink(g, i, j);
      if ( !dg_biconnected(g) )
        dg_link(g, i, j);
    } else {
      dg_link(g, i, j);
    }
    if (t >= nequil) {
      if (t == nequil) t1 = clock();
      sum += lookup ? dg_rhsc(g) : dg_rhsc_low(g);
      //sum += lookup ? dg_hsfb(g) : dg_cliquesep(g) ? 0 : dg_hsfbrjw(g);
    }
  }
  printf("star content, n %d, method %s, time used: %gs/%d\n", n,
      lookup ? "lookup" : "direct", 1.*(clock() - t1) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



int main(void)
{
  edges_t ref5[] = {
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

  cmpref(5, ref5);
  cmpref(6, ref6);
  testspeed(5, 1000000, 1);
  testspeed(6, 1000000, 1);
  testspeed(7, 1000000, 1);
  testspeed(8, 10000000, 1);
  testspeed(8, 10000, 0);
  return 0;
}
