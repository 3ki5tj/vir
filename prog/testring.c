/* iterate over permutations of n */
#include "dgring.h"
#include "dgsc.h"
#include "dgmap.h"
#include "testutil.h"



/* test the ring contents of diagrams against the reference values */
INLINE void cmpref(int n, dgref_t *ref)
{
  int i, err;
  double nr0, nr1, nr2, nr3, nr4, nr5, nr6;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr != DGREF_NPRMAX; i++) {
    dgref_build(g, ref + i);
    if ( !dg_biconnected(g) )
      continue;
    nr0 = nr1 = nr2 = nr3 = nr4 = nr5 = nr6 = 0;
    dg_fbnr_spec(g, &nr0, &err);
    if (dg_nedges(g) <= g->n + 12) {
      dgsc_iter(g, &nr1);
      dgsc_recur(g, &nr2, dg_nedges(g), 0, 0);
    } else {
      nr1 = nr2 = ref[i].nr;
      fprintf(stderr, "skipping n %d, case %d from the star content route\n", n, i);
    }
    nr3 = dgring_perm(g);
    nr4 = dgring_inv(g);
    nr5 = dgring_dp(g);
    nr6 = dgring_karp(g);
    if ( fabs(nr1 - ref[i].nr) > 1e-3 || fabs(nr2 - ref[i].nr) > 1e-3
      || fabs(nr3 - ref[i].nr) > 1e-3 || fabs(nr4 - ref[i].nr) > 1e-3
      || fabs(nr5 - ref[i].nr) > 1e-3 || fabs(nr6 - ref[i].nr) > 1e-3
      || (err == 0 && fabs(nr0 - ref[i].nr) > 1e-3) ) {
      printf("n %d: model %d nr mismatch %g,%g,%g,%g,%g,%g,%g vs %d (ref), err %d\n",
          n, i, nr0, nr1, nr2, nr3, nr4, nr5, nr6, ref[i].nr, err);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, # of subrings of %d reference diagrams verified\n", n, i);
}



static void verifyall(int n)
{
  dg_t *g;

  g = dg_open(n);
  dg_full(g);
#ifdef DGMAP_EXISTS
  if (n <= DGMAP_NMAX) {
    clock_t t0;
    int ig, err;
    double nr0, nr1, nr2;

    printf("verifying the ring contents of all diagrams\n");
    t0 = clock();
    dgmap_getuid(g); /* automatically activate the look up table */
    printf("ring content, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
    /* compare with the RJW result */
    for (ig = 0; ig < dgmap_[n].ng; ig++) {
      dgword_t code = dgmap_[n].first[ig];
      dg_decode(g, &code);
      /* skip disconnected diagrams */
      if ( !dg_biconnected(g) ) continue;
      /* check results from different methods */
      nr0 = nr1 = nr2 = 0;
      dg_fbnr_spec(g, &nr0, &err);
      nr1 = dgring_perm(g);
      nr2 = dgring_inv(g);
      if ( fabs(nr1 - nr2) > 1e-3 || (err == 0 && fabs(nr1 - nr0) > 1e-3) ) {
        printf("nr: mismatch %g, %g, %g, err %d\n", nr0, nr1, nr2, err);
        dg_print(g);
        exit(1);
      }
    }
  }
#endif
  dg_close(g);
}



static void testspeed(int n, int nsamp, int nedmax, int nedmin, int method)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, ned, eql = 1, nequil = 1000, isamp = 0;
  int goodm = 0, totm = 0, goodp = 0, totp = 0;
  double tsum = 0, sum = 0;
  double rnm = 0.9, rnp = 0.1; /* rate for moves increasing edges */

  nedmin += n * (n - 1) / 2;
  printf("speed test for n %d, nsamp %d, nedmax %d, nedmin %d\n",
      n, nsamp, nedmax, nedmin);
  dg_full(g);
  ned = dg_nedges(g);
  for (t = 1; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge0(g, &ned, nedmin, rnm, nedmax, rnp);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    if (t % 10 != 0) continue; /* avoid successive correlation */
    if (nedmax < n*(n - 1)/2)
      adjustrnp(ned, nedmax, t, 1000000, &goodp, &totp, &rnp);
    if (nedmin > 0)
      adjustrnm(ned, nedmin, t, 1000000, &goodm, &totm, &rnm);
    if ( ned > nedmax || ned < nedmin) continue;

    t0 = clock();
    if (method == 'p') {
      sum += dgring_perm(g);
    } else if (method == 'd') {
      sum += dgring_dp(g);
    } else if (method == 'i') {
      sum += dgring_inv(g);
    } else if (method == 'k') {
      sum += dgring_karp(g);
    } else {
      sum += dgring_nr(g);
    }
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("dg_ring: n %d; ave. %g; time used: %gs/%d = %gmcs\n",
      n, sum/nsamp, tsum, nsamp, tsum/nsamp*1e6);
  dg_close(g);
}



int main(int argc, char **argv)
{
  int n0 = 9, n, nsamp = 100000, nedmax = 1000000, nedmin = -1000000;
  char method = '0'; /* default */

  if (argc >= 2) n0 = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  if (argc >= 5) method = argv[4][0];
  if (argc >= 6) nedmin = atoi(argv[5]);

#ifndef N
  for (n = DGREF_NMIN; n <= DGREF_NMAX; n++) cmpref(n, dgrefs[n]);
  n = n0;
#else
  #if (N >= DGREF_NMIN && N <= DGREF_NMAX)
    cmpref(N, dgrefs[N]);
  #endif
  n = N;
#endif

  verifyall(n);

  /* with default setting, nedmax = inf, N not predefined
   * T60 timing Feb. 12 2013
   * 2.6mcs, n = 7
   * 6.6mcs, n = 8
   * 18.4mcs, n = 9
   * 53mcs, n = 10
   * 128mcs, n = 11 */
  testspeed(n, nsamp, nedmax, nedmin, method);

  DG_FREEMEMORIES()
  return 0;
}

