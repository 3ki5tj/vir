#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content by the direct method */



#include <time.h>
#include <limits.h>
#include "dg.h"
#include "dgcsep.h"



/* conversion between star-content and hard-sphere weight */
#define DG_SC2FB(sc, ned) ((sc) * (1 - (ned % 2) * 2))
#define DG_FB2SC(fb, ned) DG_SC2FB(fb, ned)



/* Ree-Hoover formula for the star content of a larger diagram of
 * n vertices from that of a smaller diagram of n-1 vertices */
static int dg_rhiter(int n, int n0, int sc)
{
  int i;

  if (n < n0) return 0;
  for (i = n0; i < n; i++) /* Ree and Hoover formula */
    sc *= (i - 1) * (i % 2 ? -1 : 1);
  return sc;
}



/* minimize the diagram of `g' by removing fully-connected vertices */
static dg_t *dg_mintop(dg_t *g)
{
  int i;
  code_t mask = mkbitsmask(g->n);

  /* remove fully-connected vertices until it is no longer biconnected */
  for (i = 0; i < g->n; i++) {
    /* construct a vertex set without i */
    code_t vs = mask ^ MKBIT(i);
    /* see if removing a fully-connected vertex leaves
     * the diagram biconnected */
    if ( g->c[i] == vs && dg_biconnectedvs(g, vs) ) {
      dg_remove1(g, g, i);
    } /* otherwise keep searching */
  }
  return g;
}



/* recursively find the star content
 * starting from the edge (i, j + 1) */
static int dg_rhsc_recur(dg_t *g, int sgn, int i, int j)
{
  int sc = 0, n = g->n;

  if (++j >= n) j = (++i) + 1;
  /* loop over pairs */
  for (; i < n - 1; j = (++i) + 1) {
    /* if the degree <= 2, removing an edge connecting i
     * makes the biconnectivity impossible */
    if (dg_deg(g, i) <= 2) continue;
    for (; j < n; j++) {
      /* try to remove the edge i-j */
      if (!dg_linked(g, i, j) || dg_deg(g, j) <= 2)
        continue;
      DG_UNLINK(g, i, j);
      if ( dg_biconnected(g) ) {
        sc += sgn /* add the diagram without (i, j) */
            + dg_rhsc_recur(g, -sgn, i, j); /* diagrams without (i, j) */
      }
      DG_LINK(g, i, j); /* link back, find subdiagrams with (i, j) */
    }
  }
  return sc;
}



#define dg_rhsc_spec(g, err) dg_rhsc_spec0(g, 0, NULL, NULL, err)

/* compute the star content (SC) of a Ree-Hoover diagram in special cases
 * only use cheap strategies to deduce the SC
 * assuming the diagram is biconnected
 * if successful, *err = 0, otherwise *err = 1
 * nocsep = 1 means the diagram has been tested with no clique separator
 *        = 0 means it MAY have clique separators
 *   if *err = 0, the diagram has a clique separator if the returned
 *                value of SC is zero
 *   if *err = 1, the diagram has no clique separator on return
 * *ned: number of edges; degs: unsorted degree sequence
 * if ned != NULL, *ned is computed on return */
INLINE int dg_rhsc_spec0(const dg_t *g, int nocsep,
    int *ned, int *degs, int *err)
{
  int i, j, n = g->n, ned0, ned1;
  static int ldegs[DG_NMAX]; /* local buffer for the degree sequence */

  *err = 0;
  /* compute the degrees of all vertices */
  if (degs == NULL) degs = ldegs;
  if (ned == NULL || *ned <= 0) {
    ned0 = dg_degs(g, degs);
    if (ned) *ned = ned0;
  } else {
    ned0 = *ned;
  }

  /* loosely connected diagrams */
  if (ned0 <= n + 1) {
    if (ned0 == n || nocsep) /* ring diagram or interwined rings (see below) */
      return 1;

    /* ned0 == n + 1
     * (a) if the two deg-3 vertices are mutually connected,
     * then the additional edge forms a clique that separates the ring
     * so the star content is zero,
     * (b) if the two deg-3 vertices are connected to deg-2 vertices
     * then no subgraph is biconnected, for removing any edge creates
     * a deg-1 vertex, making the subgraph impossible to be biconnected
     * we try to find which case by the degree sequence
     * this should be faster than computing the clique separator */
    for (i = 0; i < n; i++)
      if (degs[i] == 3) break;
    for (j = i + 1; j < n; j++)
      if (degs[j] == 3) break;
    return !dg_linked(g, i, j);
  }

  /* densely-connected diagrams */
  ned1 = n * (n - 1) /2 - ned0;
  if (ned1 < 3) {
    /* with no clique separactor, this covers the following three cases
     * `ned1' can only be 0 or 2 */
    if (nocsep) return dg_rhiter(n, ned1 + 2, 1);

    if (ned1 == 0) { /* fully connected */
      return dg_rhiter(n, 2, 1);
    } else if (ned1 == 1) { /* always has a clique separator */
      return 0;
    } else { /* ned1 == 2 */
      /* if there are two wiggly lines, the SC is nonzero only if
       * the two wiggly lines are not connected to the same vertex */
      for (i = 0; i < n; i++)
        if (degs[i] <= n - 3) return 0;
      return dg_rhiter(n, 4, 1);
    }
  }

  /* general case: try to find a clique separator */
  *err = (nocsep || !dg_cliquesep(g));
  return 0;
}



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed
 *    - # of biconnected subgraphs with odd edges removed */
static int dg_rhsc_directlow(const dg_t *g)
{
  int sc;
  dg_t *g0 = NULL;

  /* first find the minimal set of vertices that
   * contain the wiggly lines */
  g0 = dg_clone(g);
  dg_mintop(g0);
  sc = 1 + dg_rhsc_recur(g0, -1, 0, 0);
  /* use the Ree-Hoover formula to go from n0 to n */
  sc = dg_rhiter(g->n, g0->n, sc);
  dg_close(g0);
  return sc;
}



#define dg_rhsc_direct(g) dg_rhsc_direct0(g, 0, NULL, NULL)

/* directly compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed
 * nocsep = 1 means if the graph has been tested with no clique separator
 *        = 0 means it MAY have clique separators
 * *ned: number of edges; degs: degree sequence */
INLINE int dg_rhsc_direct0(const dg_t *g, int nocsep, int *ned, int *degs)
{
  int sc, err;

  /* detect special cases when possible */
  sc = dg_rhsc_spec0(g, nocsep, ned, degs, &err);
  if (err == 0) return sc;
  else return dg_rhsc_directlow(g);
}



/* compute the star content by a look up table */
INLINE int dg_rhsc_lookup(const dg_t *g)
{
  static int *sc[DGMAP_NMAX + 1]; /* biconnectivity of unique diagrams */
  int n = g->n;
  dgmap_t *m = dgmap_ + n;
  code_t c;

  if (n <= 1) return 1;
  if (sc[n] == NULL) { /* initialize the look-up table */
    dg_t *g1;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock(), t1;

    if (n >= 8) printf("n %d: initializing...\n", n);
    dgmap_init(m, n); /* compute the permutation mapping */
    if (sc[n] == NULL) xnew(sc[n], m->ng);

    t1 = clock();
    if (n >= 8) printf("n %d: diagram-map initialized %gs\n",
        n, 1.*(t1 - t0)/CLOCKS_PER_SEC);
    /* loop over unique diagrams */
    g1 = dg_open(n);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g1, &m->first[k]);
      if ( dg_biconnected(g1) ) {
        sc[n][k] = dg_rhsc_direct(g1);
        cnt++;
        nz += (sc[n][k] != 0);
      } else sc[n][k] = 0;
    }
    dg_close(g1);
    printf("n %d, computed star contents of %d/%d biconnected diagrams, %gs\n",
        n, cnt, nz, 1.*(clock() - t1)/CLOCKS_PER_SEC);
  }
  dg_encode(g, &c);
  return sc[n][ m->map[c] ]; /* m->map[c] is the id of the unique diagram */
}



#define dg_rhsc(g) dg_rhsc0(g, 0, NULL, NULL)

/* compute the star content (SC) of a Ree-Hoover diagram
 * see the comments of dg_rhsc_direct0() for details */
INLINE int dg_rhsc0(const dg_t *g, int nocsep, int *ned, int *degs)
{
  int sc, err;

  if (g->n <= DGMAP_NMAX) {
    return dg_rhsc_lookup(g);
  } else {
    sc = dg_rhsc_spec0(g, nocsep, ned, degs, &err);
    if (err == 0) return sc;
    else return dg_rhsc_directlow(g);
  }
}



#endif

