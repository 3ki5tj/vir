#include "dgcryr.h"
#include "testutil.h"



static void cmpg(dg_t *g)
{
  int i, j, k, n = g->n, err;

  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ ) {
      double fb1[3], fb2[3];
      err = 0;
      for ( k = 0; k < 3; k++ ) {
        fb1[k] = dgsc_yriter(g, i, j, k);
        fb2[k] = (double) dgrjw_yrfb(g, i, j, k);
        if ( fabs(fb1[k] - fb2[k]) > 1e-3 ) err = 1;
      }
      if ( err )
        printf("i %d, j %d, %g, %g | %g, %g | %g, %g\n",
               i, j, fb1[0], fb2[0], fb1[1], fb2[1], fb1[2], fb2[2]);
    }
}



static void test4(void)
{
  dg_t *g;

  g = dg_open(4);
  dg_full(g);
  dg_unlink(g, 0, 1);
  dg_print(g);
  cmpg(g);
  dg_close(g);
}



static void test5(void)
{
  dg_t *g;

  g = dg_open(5);
  dg_link(g, 0, 1);
  dg_link(g, 1, 2);
  dg_link(g, 2, 3);
  dg_link(g, 3, 4);
  dg_link(g, 4, 0);
  dg_link(g, 2, 4);
  dg_print(g);
  cmpg(g);
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

    cmpg(g);
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
  double fb = 0, fb1, fb2, sum = 0, tsum = 0;
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
    //if ( dg_cliquesep(g) ) continue; /* avoid clique separable */

    t0 = clock();
    for ( a = 0; a < DG_N_; a++ )
      for ( b = a + 1; b < DG_N_; b++ ) {
        if (method == 'l') {
#ifdef DGMAP_EXISTS
          fb = dgmap_yrfb(g, a, b, type);
#endif
        } else if (method == 's' || method == 'r' || method == 'i') { /* Ree-Hoover star content*/
          fb = dgsc_yriter(g, a, b, type);
        } else if (method == 'w') { /* Wheatley */
          fb = (double) dgrjw_yrfb(g, a, b, type);
        } else if (method == 'p') { /* comparison */
          fb1 = (double) dgrjw_yrfb(g, a, b, type);
          fb2 = dgsc_yriter(g, a, b, type);
          if (fabs(fb1 - fb2) > 1e-6) {
            fprintf(stderr, "fb %g (rjw) vs %g (sc)\n", fb1, fb2);
            dg_print(g);
            exit(1);
          }
          fb = fb1;
        } else { /* default */
          fb = dg_yrfb(g, a, b, type);
        }

        sum += fb;
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
  int n = 7, nsamp = 1000, type = 0, nedmax = 100000000, method = '0';

  test4();
  test5();

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) type = atoi(argv[3]);
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
