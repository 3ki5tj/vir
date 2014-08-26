#include "dgcryr.h"
#include "testutil.h"



static void cmpg(dg_t *g, int verbose)
{
  int i, j, type, n = g->n, err;
  dg_t *g1 = dg_open(n);

  dg_copy(g1, g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ ) {
      double fb0[YRTYPE_TOTAL], fb1[YRTYPE_TOTAL], fb2[YRTYPE_TOTAL];
      err = 0;
      for ( type = 0; type < 3; type++ ) {
        fb0[type] = dg_yrfb0(g, i, j, type, NULL, NULL);
        fb1[type] = dgsc_yriter(g, i, j, type);
        fb2[type] = (double) dgrjw_yrfb(g, i, j, type);
        if ( fabs(fb1[type] - fb2[type]) > 1e-3 ) err = 1;
        if ( fabs(fb0[type] - fb2[type]) > 1e-3 ) err = 2;
        die_if(dg_cmp(g1, g) != 0, "graph is corrupted, g->n %d\n", g->n);
      }
      if ( err || verbose ) {
        printf(" graph below, roots (%d, %d), "
               "type 0: %g, %g, %g | type 1: %g, %g, %g | type 2: %g, %g, %g\n",
               i, j,
               fb0[0], fb1[0], fb2[0],
               fb0[1], fb1[1], fb2[1],
               fb0[2], fb1[2], fb2[2]);
      }
      if (err) {
        dg_print(g);
        getchar();
      }
    }
  dg_close(g1);
}



static void test4(void)
{
  dg_t *g;

  g = dg_open(4);
  dg_full(g);

  dg_unlink(g, 0, 1);
  cmpg(g, 1);
  dg_print(g);

  dg_unlink(g, 2, 3);
  cmpg(g, 1);
  dg_print(g);

  dg_close(g);
}



static void test5(void)
{
  dg_t *g;

  g = dg_open(5);
  dg_link(g, 0, 2);
  dg_link(g, 0, 3);
  dg_link(g, 0, 4);
  dg_link(g, 1, 2);
  dg_link(g, 1, 4);
  dg_link(g, 2, 3);
  dg_link(g, 3, 4);
  dg_print(g);
  cmpg(g, 1);
  dg_close(g);
}



/* check if the Wheatley's method yields the same result with the direct method */
static void testconsist(int n, int nsamp)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0;
  dg_t *g;

  g = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    dg_rndswitchedge(g, &ned, 1000000, 1);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 100 != 0) continue;

    cmpg(g, 0);
    if (++isamp >= nsamp) break;
  }
  dg_close(g);
  printf("consistency check for n = %d, finished for %d samples\n", n, nsamp);
}



/* test the speed of computing fb */
static void testspeed(int n, int nsamp, int type, int nedmax, int method)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0, a, b;
  dg_t *g;
  clock_t t0;
  double fb = 0, sum = 0, tsum = 0;
  double rnp = 0.1; /* rate for moves increasing edges */

  printf("speed test for n %d, nsamp %d, nedmax %d, method %c\n",
      n, nsamp, nedmax, method);
  g = dg_open(n);
  dg_full(g);

  method = (char) (unsigned char) tolower(method);
  if (method == 'l') { /* lookup */
#ifdef DGMAP_EXISTS
    t0 = clock();
    dgmap_fb(g); /* initialization */
    fprintf(stderr, "hard-sphere weight, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
#else
    fprintf(stderr, "dgmap unavailable\n");
    exit(1);
#endif
  }

  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge(g, &ned, nedmax, rnp);

    if (eql && t >= nequil && ned <= nedmax) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 100 != 0) continue; /* avoid successive correlation */

    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if (dg_nedges(g) != ned) {
      dg_fprint(g, stderr);
      fprintf(stderr, "ned %d vs. %d\n", ned, dg_nedges(g));
      exit(1);
    }
    if ( ned > nedmax ) continue;

    t0 = clock();
    {
      for ( a = 0; a < DG_N_; a++ )
        for ( b = a + 1; b < DG_N_; b++ ) {
#ifdef DGMAP_EXISTS
          if (method == 'l') {
            fb = dgmap_yrfb(g, a, b, type);
          } else
#endif
          if (method == 's' || method == 'r' || method == 'i') { /* Ree-Hoover star content*/
            fb = dgsc_yriter(g, a, b, type);
          } else if (method == 'w') { /* Wheatley */
            fb = (double) dgrjw_yrfb(g, a, b, type);
          } else if (method == 'p') { /* comparison */
            double fb1 = (double) dgrjw_yrfb(g, a, b, type);
            double fb2 = dgsc_yriter(g, a, b, type);
            double fb3 = dgmap_yrfb(g, a, b, type);
            if (fabs(fb1 - fb3) > 1e-6 || fabs(fb1 - fb2) > 1e-3) {
              fprintf(stderr, "fb %g (rjw) vs %g (direct) vs %g (map), type %d, roots %d, %d\n",
                       fb1, fb2, fb3, type, a, b);
              dg_print(g);
              dgmap_yrfb(g, a, b, type);
              exit(1);
            }
            fb = fb1;
          } else { /* default */
            fb = dg_yrfb(g, a, b, type);
          }
          sum += fb;
        }
    }
    tsum += clock() - t0;

#ifdef VERBOSE
    printf("i %d, t %d, ned %d, fb %g\n", isamp, t, ned, fb);
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



int main(int argc, char **argv)
{
  int n = 4, nsamp = 1000, type = 0, nedmax = 100000000, method = 'p';

  test4();
  test5();

  if (argc >= 2) n = atoi(argv[1]);
  if (n < 4) n = 4;
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) type = atoi(argv[3]);
  DGYR_CHECKTYPE(type);
  if (argc >= 5) nedmax = atoi(argv[4]);
  die_if (nedmax > 0 && nedmax <= n, "nedmax %d <= n %d\n", nedmax, n);
  if (argc >= 6) method = argv[5][0];

  if ( nedmax > 0 ) {
    testspeed(n, nsamp, type, nedmax, method);
  } else {
    testconsist(n, nsamp);
  }

  DG_FREEMEMORIES()
  return 0;
}
