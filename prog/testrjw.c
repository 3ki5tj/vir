#include <time.h>
#include "dgrjw.h"



typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int fb; /* the correct weight */
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j;
  double fb;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    fb = dg_hsfb(g);
    if (!dg_biconnected(g))
      continue;
    if (fabs(fb - ref[i].fb) > 1e-3) {
      printf("n %d: model %d fb mismatch %g vs %d (ref)\n",
          n, i, fb, ref[i].fb);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, hard-sphere weights of %d reference diagrams verified\n", n, i);
}



/* test the speed of computing fb */
static void testspeed(int n, int nsamp, int nedmin, int nedmax, char method)
{
  int t, ipr, npr = n * (n - 1)/2, i, j, ned;
  double fb, sum = 0;
  int eql = 1, nequil = 1000, isamp = 0, acc;
  dg_t *g;
  clock_t t0;
  double tsum = 0;
  double rnp = 0.1; /* rate for moves increasing edges */

  printf("speed test for n %d, nsamp %d, nedmax %d, method %c\n",
      n, nsamp, nedmax, method);
  g = dg_open(n);
  dg_full(g);

  method = ctolower(method);
  if (method == 'l') { /* lookup */
    t0 = clock();
    dg_hsfb(g); /* initialization */
    printf("hard-sphere weight, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
  }

  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    /* randomly switch an edge */
    ipr = (int) (npr * rnd0());
    parsepairindex(ipr, n, &i, &j);
    if (dg_linked(g, i, j)) { /* unlink (i, j) */
      dg_unlink(g, i, j);
      acc = dg_biconnected(g);
      if ( acc ) {
        ned--;
      } else {
        dg_link(g, i, j);
      }
    } else { /* link (i, j) */
      acc = 1;
      if (ned >= nedmax) /* avoid increasing edges */
        acc = (rnd0() < rnp);
      if (acc) {
        ned++;
        dg_link(g, i, j);
      }
    }
    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */

    if (t % 1000000 == 0) { /* adjust rnp */
      rnp *= 0.5;
      printf("adjusting rnp to %g\n", rnp);
    }
    die_if (dg_nedges(g) != ned, "ned %d vs. %d\n", ned, dg_nedges(g));
    if ( ned < nedmin || ned > nedmax) continue;
    if ( dg_cliquesep(g) ) continue; /* avoid clique separable */

    t0 = clock();
    if (method == 'l') {
      /* the default function dg_hsfb() automatically invokes
       * the lookup table when possible */
      fb = dg_hsfb(g);
    } else if (method == 's' || method == 'r') { /* Ree-Hoover star content*/
      fb = dg_rhsc_direct(g) * (1 - (ned % 2) * 2);
    } else if (method == 'w') { /* Wheatley */
      fb = (double) dg_hsfb_rjw(g);
    } else { /* default */
      fb = dg_hsfb_mixed(g);
    }
    tsum += clock() - t0;
    sum += fb;
    printf("i %d, t %d, ned %d, fb %g\n", isamp, t, ned, fb);
#if 0
    if (dg_hsfb_mixed(g) != dg_hsfb_rjw(g)) {
      printf("corruption %d vs %d csep %d\n", dg_hsfb_mixed(g), dg_hsfb_rjw(g), dg_cliquesep(g));
      dg_print(g);
      exit(1);
    }
#endif
    if (++isamp >= nsamp) break;
  }

  tsum /= CLOCKS_PER_SEC;
  printf("star content, n %d, samples %d, steps %d, nedges %d-%d; "
      "method %c, fbav %g, time used: %gs/%d = %gms\n",
      n, nsamp, t, nedmin, nedmax, method,
      1.*sum/nsamp, tsum, nsamp, tsum / nsamp * 1000);
  dg_close(g);
}



int main(int argc, char **argv)
{
  int n = 12, nsamp = 1000, nedmin = 0, nedmax = 1000000;
  char method = 'd'; /* default */

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
  edges_t ref7[] = {
    {0, {{0, 0}}, -120},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {2, 3}}, 60},
    {-1, {{0, 0}}, 0},
  };

  cmpref(4, ref4);
  cmpref(5, ref5);
  cmpref(6, ref6);
  cmpref(7, ref7);

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  if (argc >= 5) method = argv[4][0];
  testspeed(n, nsamp, nedmin, nedmax, method);
  mtsave(NULL);
  return 0;
}
