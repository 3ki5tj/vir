#include "dgrjw.h"

#include <time.h>
#include "dgsc.h"


typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int fb; /* the correct weight */
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j, fb;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    fb = dg_hsfb(g);
    if (!dg_biconnected(g))
      continue;
    if (fb != ref[i].fb) {
      printf("n %d: model %d fb mismatch %d vs %d (ref)\n",
          n, i, fb, ref[i].fb);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, hard-sphere weights of %d reference diagrams verified\n", n, i);
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
    dg_hsfb(g); /* initialization */
    printf("hard-sphere weight, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
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
      /* the default function dg_hsfb() uses the lookup table
        * only if there is one */
      //if (dg_nedges(g) < n*7/3)
      sum += lookup ? dg_hsfb(g) : dg_hsfbmixed(g);
#if 0
      if (dg_hsfbmixed(g) != dg_hsfbrjw(g)) {
        printf("corruption %d vs %d csep %d\n", dg_hsfbmixed(g), dg_hsfbrjw(g), dg_cliquesep(g));
        dg_print(g);
        exit(1);
      }
#endif
    }
  }
  printf("star content, n %d, method %s, time used: %gs/%d\n", n,
      lookup ? "lookup" : "direct", 1.*(clock() - t1) / CLOCKS_PER_SEC, nsteps);
  dg_close(g);
}



int main(void)
{
  edges_t ref4[] = {
    {0, {{0, 0}}, -2},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {1, 3}}, 0},
    {2, {{0, 1}, {2, 3}}, 1},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref5[] = {
    {0, {{0, 0}}, -6},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, 3},
    {2, {{0, 1}, {1, 2}}, 0},
    {3, {{0, 1}, {2, 3}, {3, 4}}, 2},
    {3, {{0, 1}, {1, 3}, {2, 4}}, 2},
    {3, {{0, 1}, {1, 3}, {1, 4}}, 0},
    {4, {{0, 1}, {2, 3}, {3, 4}, {2, 4}}, 1},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -1},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref6[] = {
    {0, {{0, 0}}, -24},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, 12},
    {3, {{0, 1}, {2, 3}, {3, 4}}, 8},
    {4, {{0, 1}, {2, 3}, {3, 4}, {2, 4}}, 4},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -4},
    {-1, {{0, 0}}, 0},
  };

  cmpref(4, ref4);
  cmpref(5, ref5);
  cmpref(6, ref6);
  //testspeed(7, 1000000, 1);
  //testspeed(8, 10000000, 1);
  testspeed(12, 10000, 0);
  mtsave(NULL);
  return 0;
}
