#include "dgcsep.h"
#ifdef VERIFY
#include "dgmap.h"
#endif
#include "testutil.h"



/* test the RTL example */
static void testrtl(void)
{
  dg_t *g = dg_open(9), *f = dg_open(9), *f2 = dg_open(9);
  int a[9], i;
  /* RTL paper, with vertex index = 9 - i, P274 Fig. 2
   * Algorithmic aspects of vertex elimination on graphs
   * D. J. Rose, R. E. Trajan, and G. S. Lueker,
   * SIAM, J. Comput. Vol. 5. No. 2. June 1976 */
  int pairs[][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 5}, {1, 4}, {2, 3},
    {2, 6}, {3, 5}, {3, 7}, {4, 5}, {4, 8}, {5, 7}, {6, 7}, {7, 8},
    {-1, -1}};
  /* the correct fill-in graph */
  int fpairs[][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 5}, {1, 4}, {2, 3},
    {2, 6}, {3, 5}, {3, 7}, {4, 5}, {4, 8}, {5, 7}, {6, 7}, {7, 8},
    {1, 2}, {1, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 6}, {4, 6}, {4, 7}, {5, 6},
    {-1, -1}};

  dg_linkpairs(g, pairs);
  dg_print(g);
  dg_minimalorder(g, f, a);
  dg_print(f);
  dg_linkpairs(f2, fpairs);
  /* check the fill-in `f' against the correct result `f2' */
  for (i = 0; i < 9; i++) {
    if ( dgvs_neq(f->c[i], f2->c[i]) ) {
      printf("graph differs at %d\n", i);
      dgvs_printn(f->c[i], "f");
      dgvs_printn(f2->c[i], "f2");
      dg_print(f);
      dg_print(f2);
      exit(1);
    }
  }
  for (i = 0; i < 9; i++) {
    printf("%d %d\n", i, a[i]);
  }
  dg_close(g);
  dg_close(f);
  dg_close(f2);
}



static void printcsep(dg_t *g)
{
  dgvsref_t c;
  dgvs_t vs;

  c = dg_csep(g);
  if ( c ) {
    DGVS_CPY(vs, c)
    dgvs_printn(vs, NULL);
  }
  printf("\n");
}



static void testcsep(void)
{
  dg_t *g;
  /* Fig 2. of Decomposition by clique separator,
   * R. E. Tarjan, Discrete Mathematics 55 (1985) 221-232 */
  int pair3[][2] = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 4},
    {2, 4}, {3, 4}, {-1, -1}};
  /* Fig. 1 and Fig. 3 of the same paper */
  int pair4[][2] = {{0, 8}, {0, 9}, {0, 10}, {1, 2}, {1, 3}, {2, 8},
    {3, 7}, {3, 8}, {4, 5}, {4, 9}, {5, 6},  {6, 7}, {6, 9}, {8, 10},
    {9, 10}, {-1, -1}};

  g = dg_open(4);
  dg_full(g);
  dg_unlink(g, 0, 2);
  dg_print(g);
  printf("case 1: clique separator (should be 1, 3): ");
  printcsep(g);
  dg_close(g);

  g = dg_open(5);
  dg_full(g);
  dg_unlink(g, 0, 1);
  dg_unlink(g, 2, 3);
  dg_unlink(g, 3, 4);
  dg_print(g);
  printf("case 2: clique separator (should be none): ");
  printcsep(g);
  dg_close(g);

  g = dg_open(5);
  dg_linkpairs(g, pair3);
  dg_print(g);
  printf("case 3: clique separator (0, 4) or (0, 2, 4): ");
  printcsep(g);
  dg_close(g);

  g = dg_open(11);
  dg_linkpairs(g, pair4);
  dg_print(g);
  printf("case 4: clique separator (3, 8): ");
  printcsep(g);
  dg_close(g);
}



static void testfillin(void)
{
  dg_t *g, *f, *f2;
  int prg[][2] = {{0, 2}, {0, 3}, {0, 5}, {1, 3}, {1, 4}, {1, 6},
                  {3, 5}, {-1, -1}};
  int prf[][2] = {{0, 2}, {0, 3}, {0, 5}, {1, 3}, {1, 4}, {1, 6},
                  {2, 3}, {2, 5}, {3, 4}, {3, 5}, {3, 6},
                  {4, 5}, {4, 6}, {5, 6}, {-1, -1}};
  int a[7] = {0, 1, 2, 3, 4, 5, 6}, p[7] = {0, 1, 2, 3, 4, 5, 6};
  int i, j;

  g = dg_open(7);
  f = dg_open(7);
  f2 = dg_open(7);
  dg_linkpairs(g, prg);
  dg_linkpairs(f, prf);
  dg_fillin(g, f2, a, p);
  for (i = 1; i < 7; i++) {
    for (j = 0; j < i; j++) {
      if (dg_linked(f, i, j) != dg_linked(f2, i, j)) {
        dg_print(g);
        dg_print(f);
        dg_print(f2);
        fprintf(stderr, "i %d, j %d linked %d vs %d\n", i, j,
          dg_linked(f, i, j), dg_linked(f2, i, j));
      }
    }
  }
  dg_close(g);
  dg_close(f);
  dg_close(f2);
  fprintf(stderr, "passed the fill-in test\n");
}



/* test clique separators of diagrams against the reference values */
static void cmpref(int n, dgref_t *ref)
{
  int i, cs1, cs2, cs3;
  dg_t *g;

  g = dg_open(n);
  for (i = 0; ref[i].npr != DGREF_NPRMAX; i++) {
    dgref_build(g, ref + i);
    if ( !dg_biconnected(g) )
      continue;
    cs1 = (dg_cliquesep0(g, DGCSEP_LEXM) != 0);
    cs2 = (dg_cliquesep0(g, DGCSEP_MCSP) != 0);
    cs3 = (dg_cliquesep0(g, DGCSEP_MCSM) != 0);
    if ( cs1 != ref[i].cs || cs3 != ref[i].cs
      || (cs2 != cs1 && cs2 == 1) ) {
      printf("n %d, case %d, cs mismatch %d(LEX-M), %d(MCS-P), %d(MCS-M) vs %d\n",
          n, i, cs1, cs2, cs3, ref[i].cs);
      dg_print(g);
      exit(1);
    }
  }
  dg_close(g);
  printf("n %d, clique separators of %d reference diagrams verified\n", n, i);
}



static void speed_cliquesep(int n, int nsamp, int nedmax)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, k, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  double tsum[4] = {0}, sum[4] = {0};
  double rnp = 0.1; /* rate for moves increasing edges */
#define LISTSIZE 1000
  dg_t *gls[LISTSIZE] = {NULL};
  int lscnt = 0, kk;

  printf("speed test for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  for (t = 0; t < LISTSIZE; t++)
    gls[t] = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    /* randomly switch an edge */
    dg_rndswitchedge(g, &ned, nedmax, rnp);

    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      lscnt = 0;
      continue;
    }

    if ( ned <= nedmax ) {
      die_if (lscnt >= LISTSIZE, "lscnt %d, t %d\n", lscnt, t);
      dg_copy(gls[lscnt], g);
      lscnt++;
    }
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if (t % LISTSIZE != 0) continue;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[0] += (dg_csep(gls[kk]) == 0);
    tsum[0] += clock() - t0;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[1] += (dg_cliquesep0(gls[kk], 1) == 0);
    tsum[1] += clock() - t0;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[2] += (dg_cliquesep0(gls[kk], 2) == 0);
    tsum[2] += clock() - t0;

    t0 = clock();
    for (kk = 0; kk < lscnt; kk++)
      sum[3] += (dg_cliquesep0(gls[kk], 3) == 0);
    tsum[3] += clock() - t0;

    isamp += lscnt;
    lscnt = 0;
    if (isamp >= nsamp) break;
  }
  for (k = 0; k < 4; k++)
    tsum[k] /= CLOCKS_PER_SEC;
  printf("dg_cliquesep()/dg_csep(): n %d; nocsep %g%%/%g%%, %g%%\n"
      "time used: %gs/%d = %gmcs (csep), "
                 "%gs/%d = %gmcs (LEX-M), "
                 "%gs/%d = %gmcs (MCS-P), "
                 "%gs/%d = %gmcs (MCS-M)\n",
      n, 100.*sum[0]/nsamp, 100.*sum[1]/nsamp, 100.*sum[2]/nsamp,
      tsum[0], nsamp, tsum[0]/nsamp*1e6,
      tsum[1], nsamp, tsum[1]/nsamp*1e6,
      tsum[2], nsamp, tsum[2]/nsamp*1e6,
      tsum[3], nsamp, tsum[3]/nsamp*1e6);
  dg_close(g);
}


#ifdef VERIFY
#include "dgrjw.h"



#ifdef DGMAP_EXISTS
/* verify all graphs of n vertices with clique separators
 * have zero star content or fb */
static void verify_allzerofb(int n, int verbose)
{
  dg_t *g = dg_open(n);
  int nocs, cnt = 0, *visited;
  dgword_t c, npr;
  int cs1, cs2, cs3;
  double fb;
  unqid_t uid;
  int fpos = 0, fneg = 0, cntcsep = 0;

  die_if(n > DGMAP_NMAX, "n %d is too large\n", n);
  dgmap_biconnected(g); /* initialize the lookup table */
  xnew(visited, dgmap_[n].ng); /* number of unique diagrams */
  npr = (dgword_t) 1u << (n * (n - 1) / 2);
  for (c = 0; c < npr; c++) {
    dg_decode(g, &c);
    uid = dgmap_[n].map[c];
    if ( !dgmap_biconnected0(n, uid) )
      continue;
    nocs = dgmap_nocs(g);
    /* test if a graph with a clique separator necessarily
     * implies that fb is 0 */
    if (nocs == 0) { /* has at least one clique separator */
      cntcsep += 1;
      if (fabs(fb = dgmap_fb0(n, uid)) > 0.1) {
        fprintf(stderr, "n %d, fb %g, c %#x\n",
            n, fb, (unsigned) c);
        dg_print(g);
        exit(1);
      }
    } else if (verbose) {
      if (fabs(fb = dgmap_fb0(n, uid)) < 0.1) {
        if ( !visited[uid] ) {
          fprintf(stderr, "special: n %d, fb %g, nocs %d\n", n, fb, nocs);
          visited[uid] = 1;
          dg_print(g);
        }
      }
    }

    if ( dg_connected(g) ) {
      /* the maximal cardinality search MCS-M should yield
       * the same result as LEX-M (default) algorithm */
      if ( !nocs != (cs3 = (dg_csep0(g, DGCSEP_MCSM) != 0)) ) {
        dg_print(g);
        cs1 = (dg_csep0(g, DGCSEP_LEXM) != 0);
        fprintf(stderr, "csep-allzero %#x (LEX-M) vs %#x (MCS-M), nocs %d, %d\n",
          cs1, cs3, nocs, dg_ncsep(g));
        exit(1);
      }
      /* the maximal cardinality search MCS-P may yield
       * false negative results but not false positive ones */
      if ( !nocs != (cs2 = (dg_csep0(g, DGCSEP_MCSP) != 0)) ) {
        fpos += cs2;
        fneg += !cs2;
        if ( cs2 ) {
          dg_print(g);
          fprintf(stderr, "csep %#x (LEX-M) vs %#x (MCS-P), ncsep %d, %d\n",
            dg_csep0(g, 1), cs2, nocs, dg_ncsep(g));
          exit(1);
        }
      }
    }
    cnt++;
  }
  dg_close(g);
  printf("verified all zero fb for %d graphs of n = %d\n", cnt, n);
  printf("MCS-P: n %d, false positive %d, false negative %d, %d/%d\n",
      n, fpos, fneg, cntcsep, cnt);
}
#endif /* defined(DGMAP_EXISTS) */



/* verify graphs with clique separators have zero star content or fb */
static void verify_zerofb(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  int t, i, j, cnt = 0, fpos = 0, fneg = 0;
  int cs1, cs2, cs3;
  double fb;

  dg_full(g);
  for (t = 1; t <= nsteps * 10; t++) {
    i = randpair(n, &j);
    if (dg_linked(g, i, j)) {
      dg_unlink(g, i, j);
      if (!dg_biconnected(g)) {
        dg_link(g, i, j);
      }
    } else {
      dg_link(g, i, j);
    }
    if (t % 10 != 0) continue;
    /* test if a graph with a clique separator necessarily
     * implies that fb is 0 */
    if ( (cs1 = (dg_csep(g) != 0)) ) {
      if (fabs(fb = dg_fb(g)) > 0.1) {
        fprintf(stderr, "n %d, fb %g, cs1 %#x\n",
            n, fb, (unsigned) cs1);
        dg_print(g);
        exit(1);
      }
    }
    cs2 = (dg_csep0(g, DGCSEP_MCSP) != 0);
    if ((cs2 != 0) != (cs1 != 0)) {
      fpos += (cs2 != 0);
      fneg += (cs2 == 0);
      if (cs2 != 0) { /* false positive */
        dg_print(g);
        fprintf(stderr, "csep %#x (LEX-M) vs %#x (MCS-P) \n",
            dg_csep0(g, 1), cs2);
      }
    }
    cs3 = (dg_csep0(g, DGCSEP_MCSM) != 0);
    if ((cs3 != 0) != (cs1 != 0)) {
      dg_print(g);
      fprintf(stderr, "csep %#x (LEX-M) vs %#x (MCS-M)\n",
        dg_csep0(g, 1), cs3);
      exit(1);
    }
    cnt++;
  }
  dg_close(g);
  printf("verified zero fb for %d graphs of n = %d, LEX-P fneg %g%% (%d)\n",
      cnt, n, 100.*fneg/cnt, fneg);
}
#endif /* defined(VERIFY) */



int main(int argc, char **argv)
{
  int n = 9, nsteps = 10000000, nedmax = 100000000;

#ifdef N
  n = N;
#else
  testrtl();
  testcsep();
  testfillin();
  for (n = DGREF_NMIN; n <= DGREF_NMAX; n++) cmpref(n, dgrefs[n]);
#endif

  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsteps = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  /* timing of dg_cliquesep() LEX-M on T60 with the default setting (Dec. 31, 2013)
   * n   generic    with -DN=n
   * 9    1.7mcs      1.6mcs
   * 12   2.9mcs      2.8mcs
   * 16   5.0mcs      4.7mcs
   * */
  speed_cliquesep(n, nsteps, nedmax);

#ifdef VERIFY
  printf("\n\n\n");

  #ifdef N
  if (N <= DGMAP_NMAX)
    verify_allzerofb(N, 0);
  #else
  {
    int i;
    for (i = 3; i <= DGMAP_NMAX; i++)
      verify_allzerofb(i, 0);
  }
  #endif /* defined(N) */

  if (n > DGMAP_NMAX)
    verify_zerofb(n, nsteps);
#endif /* VERIFY */

  return 0;
}

