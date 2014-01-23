#include "dgmap.h"
#include "testutil.h"



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, dgref_t *ref)
{
  int i;
  double fb, nr;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr != DGREF_NPRMAX; i++) {
    dgref_build(g, ref + i); /* build the reference graph */
    if ( !dg_biconnected(g) )
      continue;
    nr = -1;
    fb = (double) dg_fbnr(g, &nr);
    if (fabs(fb - ref[i].fb) > 1e-3 || fabs(nr - ref[i].nr) > 1e-3) {
      printf("n %d: model %d mismatch fb %g vs %d (ref), nr %g vs %d\n",
          n, i, fb, ref[i].fb, nr, ref[i].nr);
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
  double fb, fb1, fb2, sum = 0, tsum = 0;
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

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */

    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if (dg_nedges(g) != ned) {
      dg_fprint(g, stderr);
      fprintf(stderr, "ned %d vs. %d\n", ned, dg_nedges(g));
      exit(1);
    }
    if ( ned > nedmax ) continue;
    if ( dg_cliquesep(g) ) continue; /* avoid clique separable */

    t0 = clock();
    if (method == 'l') {
#ifdef DGMAP_EXISTS
      fb = dgmap_fb(g);
#endif
    } else if (method == 's' || method == 'r' || method == 'i') { /* Ree-Hoover star content*/
      fb = dgsc_fb0(g, (method == 'r' ? DGSC_RECUR : DGSC_ITER), NULL, NULL);
    } else if (method == 'w') { /* Wheatley */
      fb = (double) dgrjw_fb(g);
    } else if (method == 'p') { /* comparison */
      fb1 = (double) dgrjw_fb(g);
      fb2 = dgsc_fb(g);
      if (fabs(fb1 - fb2) > 1e-6) {
        fprintf(stderr, "fb %g (rjw) vs %g (sc)\n", fb1, fb2);
        dg_print(g);
        exit(1);
      }
      fb = fb1;
    } else { /* default */
      fb = dg_fb(g);
    }
    tsum += clock() - t0;
    sum += fb;
#ifdef VERBOSE
    printf("i %d, t %d, ned %d, fb %g\n", isamp, t, ned, fb);
#endif
#if 0
    if (dg_fb(g) != dgrjw_fb(g)) {
      printf("corruption %d vs %d csep %d\n", dg_fb(g), dgrjw_fb(g), dg_cliquesep(g));
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



/* test the accuray of computing fb
 * compile the program using different flags -DRJW32, -DRJW64, and -DRJWDBL
 * and compare the resulting fbrjwXX.dat */
static void testaccuracy(int n, int nsamp, int nedmin)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0, err;
  dg_t *g;
  clock_t t0;
  double fb0, fb, sum = 0, tsum = 0;
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
    fb = (double) dgrjw_fb(g);
    tsum += clock() - t0;

    fb0 = dg_fbnr_spec(g, NULL, &err);
    if (err == 0) {
      die_if (fabs(fb - fb0) > 0.01, "fb %g, fb1 %g\n", fb, fb0);
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
  int n = 12, nsamp = 1000, nedmax = 1000000;
  char method = 'd'; /* default */

#ifdef N
  n = N;
#else
  int i;
#endif

#ifndef N
  for (i = DGREF_NMIN; i <= DGREF_NMAX; i++)
    cmpref(i, dgrefs[i]);
#elif (N >= DGREF_NMIN) && (N <= DGREF_NMAX)
  cmpref(N, dgrefs[N]);
#endif

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  die_if (nedmax > 0 && nedmax <= n, "nedmax %d <= n %d\n", nedmax, n);
  if (argc >= 5) method = argv[4][0];

  if (nedmax > 0) { /* do the speed test */
    /* 11.3ms by default setting T60 (Jan. 8, 2014) */
    testspeed(n, nsamp, nedmax, method);
    //mtsave(NULL);
  } else { /* do the accuracy test for Wheatley's method */
    int nedmin = n * (n - 1) / 2 + nedmax;
    testaccuracy(n, nsamp, nedmin);
  }

  dgrjw_free();
  return 0;
}
