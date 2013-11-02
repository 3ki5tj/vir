#include "dgcsep.h"
#include "testutil.h"


/* test the RTL example */
static void testrtl(void)
{
  dg_t *g = dg_open(9), *f = dg_open(9);
  int a[9], i;
  /* RTL paper index = 9 - i, P274 Fig. 2
   * Algorithmic aspects of vertex elimination on graphs
   * D. J. Rose, R. E. Trajan, and G. S. Lueker,
   * SIAM, J. Comput. Vol. 5. No. 2. June 1976 */
  int pair[][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 5}, {1, 4}, {2, 3}, {2, 6},
    {3, 5}, {3, 7}, {4, 5}, {4, 8}, {5, 7}, {6, 7}, {7, 8}, {-1, -1}};

  dg_linkpairs(g, pair);
  dg_print(g);
  dg_minimalorder(g, f, a);
  dg_print(f);
  for (i = 0; i < 9; i++) {
    printf("%d %d\n", i, a[i]);
  }
  dg_close(g);
  dg_close(f);
}



static void testcsep(void)
{
  dg_t *g;
  dgword_t c, bw;
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
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n\n");
  dg_close(g);

  g = dg_open(5);
  dg_full(g);
  dg_unlink(g, 0, 1);
  dg_unlink(g, 2, 3);
  dg_unlink(g, 3, 4);
  dg_print(g);
  printf("case 2: clique separator (should be none): ");
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n\n");
  dg_close(g);

  g = dg_open(5);
  dg_linkpairs(g, pair3);
  dg_print(g);
  printf("case 3: clique separator (0, 4) or (0, 2, 4): ");
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n\n");
  dg_close(g);

  g = dg_open(11);
  dg_linkpairs(g, pair4);
  dg_print(g);
  printf("case 4: clique separator: ");
  for (c = dg_cliquesep(g); c; c ^= bw)
    printf("%d ", bitfirstlow(c, &bw));
  printf("\n\n");
  dg_close(g);
}



static void speed_cliquesep(int n, int nsamp, int nedmax)
{
  dg_t *g = dg_open(n);
  clock_t t0;
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  double tsum = 0, sum = 0;
  double rnp = 0.1; /* rate for moves increasing edges */

  printf("speed test for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  dg_full(g);
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
    if ( ned > nedmax ) continue;

    t0 = clock();
    sum += (dg_cliquesep(g) == 0);
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("dg_cliquesep: n %d; sum %g; time used: %gs/%d = %gmcs\n",
      n, sum/nsamp, tsum, nsamp, tsum/nsamp*1e6);
  dg_close(g);
}


#ifdef VERIFY
#include "dgrjw.h"

/* verify all graphs of n vertices with clique separators
 * have zero star content or fb */
static void verify_allzerofb(int n, int verbose)
{
  dg_t *g = dg_open(n);
  int ncsep, cnt = 0, *visited;
  dgword_t c, npr;
  double fb;
  unqid_t uid;

  die_if(n > DGMAP_NMAX, "n %d is too large\n", n);
  dg_biconnected_lookup(g); /* initialize the lookup table */
  xnew(visited, dgmap_[n].ng); /* number of unique diagrams */
  npr = (dgword_t) 1u << (n * (n - 1) / 2);
  for (c = 0; c < npr; c++) {
    dg_decode(g, &c);
    uid = dgmap_[n].map[c];
    if ( !dg_biconnected_lookuplow(n, uid) )
      continue;
    ncsep = dg_ncsep_lookuplow(g, c);
    /* test if a graph with a clique separator necessarily
     * implies that fb is 0 */
    if (ncsep > 0) { /* has at least one clique separator */
      if (fabs(fb = dg_hsfb_lookuplow(n, uid)) > 0.1) {
        fprintf(stderr, "n %d, fb %g, c %#x\n",
            n, fb, (unsigned) c);
        dg_print(g);
        exit(1);
      }
    } else if (verbose) {
      if (fabs(fb = dg_hsfb_lookuplow(n, uid)) < 0.1) {
        if ( !visited[uid] ) {
          fprintf(stderr, "special: n %d, fb %g, ncsep %d\n", n, fb, ncsep);
          visited[uid] = 1;
          dg_print(g);
        }
      }
    }
    cnt++;
  }
  dg_close(g);
  printf("verified all zero fb for %d graphs of n = %d\n", cnt, n);
}




/* verify graphs with clique separators have zero star content or fb */
static void verify_zerofb(int n, int nsteps)
{
  dg_t *g = dg_open(n);
  int t, i, j, cnt = 0;
  dgword_t c;
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
    if ((c = dg_cliquesep(g)) != 0) {
      if (fabs(fb = dg_hsfb_mixed(g)) > 0.1) {
        fprintf(stderr, "n %d, fb %g, c %#x\n",
            n, fb, (unsigned) c);
        dg_print(g);
        exit(1);
      }
    }
    cnt++;
  }
  dg_close(g);
  printf("verified zero fb for %d graphs of n = %d\n", cnt, n);
}
#endif /* defined(VERIFY) */



int main(int argc, char **argv)
{
  int n = 9, nsteps = 1000000, nedmax = 1000;
#ifdef N
  n = N;
#endif
#ifndef N
  testrtl();
  testcsep();
#endif
  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsteps = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  /* T60 3mcs default setting */
  speed_cliquesep(n, nsteps, nedmax);

#ifdef VERIFY
  printf("\n\n\n");
#ifdef N
  verify_allzerofb(N, 0);
#else
  {
    int i;
    for (i = 3; i <= DGMAP_NMAX; i++)
      verify_allzerofb(i, 0);
  }
#endif
  verify_zerofb(n, nsteps);
#endif /* VERIFY */
  return 0;
}

