#define DGMAP_NEEDSSC 1
#include "dgmap.h"
#include "dgsc.h"
#include "dgrjw.h"
#include "testutil.h"



/* test the star contents of diagrams against the reference values */
static void cmpref(int n, dgref_t *ref)
{
  int i, err;
  double fb0, nr0, fb1, nr1, fb2, nr2, fb3, nr3;
  dg_t *g;

  g = dg_open(n);
  printf("n %d: ", n);
  for (i = 0; ref[i].npr != DGREF_NPRMAX; i++) {
    printf("%d ", i);
    dgref_build(g, ref + i); /* build the reference graph */
    if (!dg_biconnected(g)) {
      fprintf(stderr, "test case %d of n = %d is not connected\n", i, n);
      continue;
    }
    nr0 = nr1 = nr2 = -1;
    fb0 = dg_fbnr_spec0(g, &nr0, NULL, NULL, &err);
    fb1 = dgsc_fbnr0(g, &nr1, DGSC_ITER, NULL, NULL);
    fb2 = dgsc_fbnr0(g, &nr2, DGSC_RECUR, NULL, NULL);
    fb3 = (double) dgrjw_fb(g);
    nr3 = dgring_nr(g);
    if ( (err == 0 && (fabs(fb0 - ref[i].fb) > 0.001 || fabs(nr0 - ref[i].nr) > 0.001))
      || fabs(fb1 - ref[i].fb) > 0.001 || fabs(nr1 - ref[i].nr) > 0.001
      || fabs(fb2 - ref[i].fb) > 0.001 || fabs(nr2 - ref[i].nr) > 0.001
      || fabs(fb3 - ref[i].fb) > 0.001 || fabs(nr3 - ref[i].nr) > 0.001) {
      fprintf(stderr, "n %d: model %d mismatch, "
          "fb: %g, %g, %g, %g vs %d (ref), nr: %g, %g, %g, %g vs %d (ref)\n",
          n, i, fb0, fb1, fb2, fb3, ref[i].fb, nr0, nr1, nr2, nr3, ref[i].nr);
      dg_fprint(g, stderr);
      exit(1);
    }
  }
  dg_close(g);
  printf("\nn %d, star contents of %d reference diagrams verified\n", n, i);
}



static void verifyall(int n)
{
  dg_t *g = dg_open(n);
  dg_full(g);

#if DGMAP_EXISTS
  if (n <= DGMAP_NMAX) {
    int ig, err;
    double fb0, fb1, fb2, fb3, nr0, nr1, nr2, nr3;
    clock_t t0;

    printf("verifying the star contents of all diagrams\n");
    t0 = clock();
    dgmap_getuid(g); /* automatically activate the look up table */
    printf("star content, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
    /* compare with the RJW result */
    g = dg_open(g->n);
    for (ig = 0; ig < dgmap_[n].ng; ig++) {
      dgword_t code = dgmap_[n].first[ig];
      dg_decode(g, &code);
      /* skip disconnected diagrams */
      if ( !dg_connected(g) ) continue;
      /* check results from different methods */
      fb0 = dg_fbnr_spec(g, &nr0, &err);
      fb1 = dgsc_fbnr0(g, &nr1, DGSC_ITER, NULL, NULL);
      fb2 = dgsc_fbnr0(g, &nr2, DGSC_RECUR, NULL, NULL);
      fb3 = (double) dgrjw_fb(g);
      nr3 = dgring_nr(g);
      if ( fabs(fb1 - fb2) > 1e-3 || fabs(fb1 - fb3) > 1e-3
          || (err == 0 && fabs(fb1 - fb0) > 1e-3) ) {
        printf("fb: mismatch %g, %g, %g, %g, err %d\n", fb0, fb1, fb2, fb3, err);
        dg_print(g);
        exit(1);
      }
      if ( fabs(nr1 - nr2) > 1e-3 || fabs(nr1 - nr3) > 1e-3
          || (err == 0 && fabs(nr1 - nr0) > 1e-3) ) {
        printf("nr: mismatch %g, %g, %g, %g, err %d\n", nr0, nr1, nr2, nr3, err);
        dg_print(g);
        exit(1);
      }
    }
  }
#endif
  dg_close(g);
}

/* test the speed of computing star content */
static void testspeed(int n, int nsamp, int nedmax, char method)
{
  int t, ned, nequil = 1000, eql = 1, isamp = 0, good = 0, tot = 0;
  dg_t *g;
  clock_t t0;
  double sum = 0, tsum = 0, rnp = 0.1;

  g = dg_open(n);
  dg_full(g);
#if DGMAP_EXISTS
  if (method == 'l' && n < DGMAP_NMAX) {
    int ig, err;
    double fb0, fb1, fb2, fb3, nr0, nr1, nr2, nr3;
    dg_t *g2;

    printf("verifying the ring contents of all diagrams\n");
    t0 = clock();
    dgmap_getuid(g); /* automatically activate the look up table */
    printf("star content, n %d, initialization: %gs\n",
      n, 1.*(clock() - t0) / CLOCKS_PER_SEC);
    /* compare with the RJW result */
    g2 = dg_open(g->n);
    for (ig = 0; ig < dgmap_[n].ng; ig++) {
      dgword_t code = dgmap_[n].first[ig];
      dg_decode(g2, &code);
      /* skip disconnected diagrams */
      if ( !dg_connected(g2) ) continue;
      /* check results from different methods */
      fb0 = dg_fbnr_spec(g2, &nr0, &err);
      fb1 = dgsc_fbnr0(g2, &nr1, DGSC_ITER, NULL, NULL);
      fb2 = dgsc_fbnr0(g2, &nr2, DGSC_RECUR, NULL, NULL);
      fb3 = (double) dgrjw_fb(g2);
      nr3 = dgring_nr(g2);
      if ( fabs(fb1 - fb2) > 1e-3 || fabs(fb1 - fb3) > 1e-3
          || (err == 0 && fabs(fb1 - fb0) > 1e-3) ) {
        printf("fb: mismatch %g, %g, %g, %g, err %d\n", fb0, fb1, fb2, fb3, err);
        dg_print(g2);
        exit(1);
      }
      if ( fabs(nr1 - nr2) > 1e-3 || fabs(nr1 - nr3) > 1e-3
          || (err == 0 && fabs(nr1 - nr0) > 1e-3) ) {
        printf("nr: mismatch %g, %g, %g, %g, err %d\n", nr0, nr1, nr2, nr3, err);
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
    sum += dgsc_fbnr0(g, NULL, (method == 'r' ? DGSC_RECUR : DGSC_ITER), NULL, NULL);
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
  int i, n, nedmax = 10000000;

#ifndef N
  /* n = 8 cases are difficult for direct methods, so we stop at n = 7 */
  for (i = DGREF_NMIN; i <= DGREF_NMAX; i++)
    cmpref(i, dgrefs[i]);
  for (i = 4; i <= DGMAP_NMAX; i++)
    verifyall(i);
  n = 8;
#else
#if N <= DGREF_NMAX
  if (dgrefs[N] != NULL) cmpref(N, dgrefs[N]);
#endif
  n = N;
#endif
  verifyall(n);
  testspeed(n, 100000,  nedmax, 'd');
  DG_FREEMEMORIES()
  return 0;
}
