#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content
 * deprecated, use dgrjw.h instead */



#include <time.h>
#include <limits.h>
#include "dg.h"
#include "dgcsep.h"



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

  /* remove fully-connected vertices until it is no longer biconnected */
  for (i = 0; i < g->n; i++) {
    /* construct a vertex set without i */
    code_t vs = (((code_t) 1u << g->n) - 1) ^ ((code_t) 1u << i);
    /* see if removing a fully-connected vertex leaves
     * the diagram biconnected */
    if ( g->c[i] == vs && dg_biconnectedvs(g, vs) ) {
      dg_shrink1(g, g, i);
      g->n--;
    } /* otherwise keep searching */
  }
  return g;
}



/* recursively find the star content */
static int dg_rhsc_recur(dg_t *g, int par, int i, int j)
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
      dg_unlink(g, i, j);
      if ( dg_biconnected(g) ) {
        sc += !par ? -1 : 1; /* add the diagram without (i, j) */
        sc += dg_rhsc_recur(g, !par, i, j); /* diagrams without (i, j) */
      }
      dg_link(g, i, j); /* link back, find subdiagrams with (i, j) */
    }
  }
  return sc;
}



#define dg_rhsc_low(g) dg_rhsc_low0(g, 0)

/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed
 * nocsep: if the graph has been tested with no clique separator */
static int dg_rhsc_low0(const dg_t *g, int nocsep)
{
  int par = 0, sc, i, n = g->n, ned, ned0;
  dg_t *g0 = NULL;

  ned0 = dg_nedges(g);
  ned = n * (n - 1) /2 - ned0;
  if (ned >= 3) { /* general case */
    /* try to find a clique separator */
    if (!nocsep && dg_cliquesep(g)) return 0;

    /* if ned0 == n, it's the ring diagram
     * if ned0 == n + 1, and if the two deg-3 vertices are mutually connected,
     * then the additional edge forms a clique that separates the ring
     * and this case has been ruled out
     * or if the two deg-3 vertices are connected to deg-2 vertices
     * then no subgraph is biconnected, for removing any edge creates
     * a deg-1 vertex, making the subgraph impossible to be biconnected */
    if (ned0 <= n + 1) return 1;
    /* general case, first find the minimal set of vertices that
     * contain the wiggly lines */
    g0 = dg_clone(g);
    dg_mintop(g0);
    sc = (par ? -1 : 1) + dg_rhsc_recur(g0, par, 0, 0);
    /* use the Ree-Hoover formula to go from n0 to n */
    sc = dg_rhiter(n, g0->n, sc);
    dg_close(g0);
    return sc;
  } else if (nocsep) {
    /* with no clique separactor, this covers the following three cases */
    return dg_rhiter(n, ned + 2, 1);
  } else if (ned == 0) {
    return dg_rhiter(n, 2, 1);
  } else if (ned == 1) {
    return 0;
  } else { /* ned == 2 */
    /* if there are two wiggly lines, the SC is nonzero only if
     * the two wiggly lines are not connected to the same vertex */
    for (i = 0; i < n; i++)
      if (dg_deg(g, i) <= n - 3) return 0;
    return dg_rhiter(n, 4, 1);
  }
}



/* compute the star content by a look up table */
INLINE int dg_rhsc_lookup(const dg_t *g)
{
  static int *sc[DGMAP_NMAX + 1]; /* biconnectivity of unique diagrams */
  int n = g->n;
  dgmap_t *m = dgmap_ + n;
  code_t c;

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
        sc[n][k] = dg_rhsc_low(g1);
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



#define dg_rhsc(g) dg_rhsc0(g, 0)

/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed
 * nocsep: if the graph has been tested with no clique separator */
INLINE int dg_rhsc0(const dg_t *g, int nocsep)
{
  return (g->n <= DGMAP_NMAX) ? dg_rhsc_lookup(g) : dg_rhsc_low0(g, nocsep);
}



#endif
