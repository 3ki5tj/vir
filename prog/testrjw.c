#include "dgrjw.h"
#include "testutil.h"


typedef struct {
  int npr; /* number of wiggly lines */
  int id[64][2]; /* vertex pairs in wiggly lines */
  int fb; /* the correct weight */
} edges_t;



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, edges_t *ref)
{
  int i, j;
  double fb, fb1, fb2;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr >= 0; i++) {
    dg_full(g);
    for (j = 0; j < ref[i].npr; j++)
      dg_unlink(g, ref[i].id[j][0], ref[i].id[j][1]);
    fb = dg_hsfb(g);
    if ( !dg_biconnected(g) )
      continue;
    if (fabs(fb - ref[i].fb) > 1e-3) {
      printf("n %d: model %d fb mismatch %g vs %d (ref)\n",
          n, i, fb, ref[i].fb);
      dg_print(g);
      exit(1);
    }
    fb1 = 1.*dg_hsfb_rjw(g);
    fb2 = DG_SC2FB(dg_rhsc_directlow(g), dg_nedges(g));
    if (fabs(fb1 - fb) > 1e-3 || fabs(fb2 - fb) > 1e-3) {
      printf("n %d: model %d fb mismatch %g(rjw), %g(sc) vs %g (ref)\n",
          n, i, fb1, fb2, fb);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, hard-sphere weights of %d reference diagrams verified\n", n, i);
}



/* test the speed of computing fb */
static void testspeed(int n, int nsamp, int nedmax, char method)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  dg_t *g;
  clock_t t0;
  double fb, sum = 0, tsum = 0;
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
    dg_rndswitchedge(g, &ned, nedmax, rnp);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */

    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    die_if (dg_nedges(g) != ned, "ned %d vs. %d\n", ned, dg_nedges(g));
    if ( ned > nedmax ) continue;
    if ( dg_cliquesep(g) ) continue; /* avoid clique separable */

    t0 = clock();
    if (method == 'l') {
      /* the default function dg_hsfb() automatically invokes
       * the lookup table when possible */
      fb = dg_hsfb(g);
    } else if (method == 's' || method == 'r') { /* Ree-Hoover star content*/
      fb = DG_SC2FB( dg_rhsc_directlow(g), ned );
    } else if (method == 'w') { /* Wheatley */
      fb = (double) dg_hsfb_rjw(g);
    } else if (method == 'p') { /* comparison */
      double fb1 = (double) dg_hsfb_rjw(g);
      double fb2 = DG_SC2FB( dg_rhsc_directlow(g), ned );
      if (fabs(fb1 - fb2) > 1e-6) {
        fprintf(stderr, "fb %g (rjw) vs %g (sc)\n", fb1, fb2);
        dg_print(g);
        exit(1);
      }
      fb = fb1;
    } else { /* default */
      fb = dg_hsfb_mixed(g);
    }
    tsum += clock() - t0;
    sum += fb;
#ifdef VERBOSE
    printf("i %d, t %d, ned %d, fb %g\n", isamp, t, ned, fb);
#endif
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
  printf("star content, n %d, samples %d, steps %d, nedmax %d; "
      "method %c, fbav %g, time used: %gs/%d = %gms\n",
      n, nsamp, t, nedmax, method,
      1.*sum/nsamp, tsum, nsamp, tsum / nsamp * 1000);
  dg_close(g);
}



/* test the accuray of computing fb */
static void testaccuracy(int n, int nsamp, int nedmin)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0, err;
  dg_t *g;
  clock_t t0;
  double fb, sc, sum = 0, tsum = 0;
  double rnm = 0.1; /* rate for moves increasing edges */
  char fn[64] = "fbrjwdbl.dat";
  FILE *fp;

#ifndef DGRJW_DOUBLE
  sprintf(fn, "fbrjw%d.dat", (int) sizeof(dgrjw_fb_t) * 8);
#endif
  xfopen(fp, fn, "w", return);
  printf("accuracy test for n %d, nsamp %d, nedmin %d, size(fb_t) %d\n",
      n, nsamp, nedmin, (int) sizeof(dgrjw_fb_t));
  g = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge0(g, &ned, nedmin, rnm, 1000000, 1);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */

    adjustrnm(ned, nedmin, t, 1000000, &good, &tot, &rnm);
    die_if (dg_nedges(g) != ned, "ned %d vs. %d\n", ned, dg_nedges(g));
    if ( ned < nedmin ) continue;
    if ( dg_cliquesep(g) ) continue; /* avoid clique separable */

    t0 = clock();
    fb = (double) dg_hsfb_rjw(g);
    tsum += clock() - t0;

    sc = dg_rhsc_spec(g, &err);
    if (err == 0) {
      double fb1 = DG_SC2FB(sc, ned);
      die_if (fabs(fb - fb1) > 0.01, "fb %g, fb1 %g\n", fb, fb1);
    }
    fprintf(fp, "%.0f %d\n", fb, ned);
    if (++isamp >= nsamp) break;
  }

  tsum /= CLOCKS_PER_SEC;
  printf("Wheatley's method: n %d, samples %d, steps %d, nedmin %d; "
      "output %s, fbav %g, time used: %gs/%d = %gms\n",
      n, nsamp, t, nedmin, fn,
      1.*sum/nsamp, tsum, nsamp, tsum / nsamp * 1000);
  printf("try to compile the code with different options -DRJW32/-DRJW64/-DRJWDBL\n");
  dg_close(g);
  fclose(fp);
}



int main(int argc, char **argv)
{
  int i, n = 12, nsamp = 1000, nedmax = 1000000;
  char method = 'd'; /* default */
#ifdef N
  n = N;
#endif
  edges_t ref4[] = {
    {0, {{0, 0}}, -2},
    {1, {{0, 1}}, 0},
    {2, {{0, 1}, {1, 3}}, 0},
    {2, {{0, 1}, {2, 3}}, 1},
    {-1, {{0, 0}}, 0},
  };
  edges_t ref5[] = {
    {1, {{0, 1}}, 0},
    {0, {{0, 0}}, -6},
    {2, {{0, 1}, {2, 3}}, 3},
    {2, {{0, 1}, {1, 2}}, 0},
    {3, {{0, 1}, {0, 3}, {2, 4}}, 2},
    {3, {{0, 1}, {2, 3}, {3, 4}}, 2},
    {3, {{0, 1}, {1, 3}, {2, 4}}, 2},
    {3, {{0, 1}, {1, 3}, {1, 4}}, 0},
    {4, {{0, 1}, {2, 3}, {3, 4}, {2, 4}}, 1},
    {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -1},
    {5, {{0, 3}, {0, 4}, {1, 2}, {1, 4}, {2, 3}}, -1},
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
  edges_t *refs[8] = {NULL, NULL, NULL, NULL, ref4, ref5, ref6, ref7};

#ifndef N
  for (i = 4; i <= 7; i++)
    cmpref(i, refs[i]);
#elif (N >= 4) && (N <= 7)
  cmpref(N, refs[N]);
#endif

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  die_if (nedmax > 0 && nedmax <= n, "nedmax %d <= n %d\n", nedmax, n);
  if (argc >= 5) method = argv[4][0];

  if (nedmax > 0) { /* do the speed test */
    /* 18ms by default setting */
    testspeed(n, nsamp, nedmax, method);
    mtsave(NULL);
  } else { /* do the accuracy test for Wheatley's method */
    int nedmin = n * (n - 1) / 2 + nedmax;
    testaccuracy(n, nsamp, nedmin);
  }
  return 0;
}
