#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content */



#include <time.h>
#include <limits.h>
#include "dg.h"



/* Ree-Hoover formula for the star content of a larger diagram of
 * n vertices from that of a smaller diagram of n-1 vertices */ 
static int dg_iter(int n, int n0, int sc)
{
  int i;

  if (n < n0) return 0;
  for (i = n0; i < n; i++) /* Ree and Hoover formula */
    sc *= (i - 1) * (i % 2 ? -1 : 1);
  return sc;
}



/* minimal diagram with the same topology as `g', fully-connected
 * vertices, those connect all other vertices are removed */
static dg_t *dg_mintop(const dg_t *g)
{
  int i;
  dg_t *g0, *g1;

  g0 = dg_clone(g);
  g1 = dg_open(g->n - 1);
  /* remove vertices until it is no longer biconnected */
  while (1) {
    g1->n = g0->n - 1;
    for (i = 0; i < g0->n; i++)
      /* see if removing a fully-connected vertex leaves the diagram biconnected */
      if ( dg_deg(g0, i) == g0->n - 1
        && dg_biconnected( dg_shrink1(g1, g0, i) ) ) {
        g0->n = g1->n;
        dg_copy(g0, g1);
        break;
      } /* otherwise keep searching */
    if (i >= g0->n) break;
  }
  dg_close(g1);
  return g0;
}



/* recursively find the star content
 * edges in `c' are not to be removed */
static int dg_rhsc_recur(dg_t *a, dg_t *c, int par, int lookup)
{
  int i, j, sc, n = a->n, bc;

  for (i = 0; i < n - 1; i++) {
    /* if the degree <= 2, removing the an edge connecting i
     * makes the diagram impossible to biconnected */
    if (dg_deg(a, i) <= 2) continue;
    for (j = i + 1; j < n; j++) {
      /* try to remove the edge i-j */
      if (!dg_linked(a, i, j) || dg_linked(c, i, j) || dg_deg(a, j) <= 2)
        continue;
      dg_unlink(a, i, j);
      bc = lookup ? dg_biconnected_lookup(a) : dg_biconnected(a);
      if ( !bc ) {
        dg_link(a, i, j);
        continue;
      }
     
      /* found an edge */
      sc = !par ? -1 : 1; /* add the diagram without (i, j) */
      sc += dg_rhsc_recur(a, c, !par, lookup); /* diagrams without (i, j) */
      dg_link(a, i, j); /* link back (i, j) */

      dg_link(c, i, j); /* try diagrams must have (i, j); */
      sc += dg_rhsc_recur(a, c, par, lookup);
      dg_unlink(c, i, j); /* release the restriction */
      return sc;
    }
  }
  return 0;
}



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed */
static int dg_rhsc_low(dg_t *g, int lookup)
{
  int par = 0, sc, i, n = g->n, ned;
  dg_t *g0 = NULL, *gc = NULL;

  ned = n * (n - 1) /2 - dg_nedges(g);
  if (ned == 0) return dg_iter(n, 2, 1);
  else if (ned == 1) return 0;
  else if (ned == 2) {
    /* if there are two wiggly lines, the SC is nonzero only if
     * the two wiggly lines are not connected to the same vertex */
    for (i = 0; i < n; i++)
      if (dg_deg(g, i) <= n - 3) return 0;
    return dg_iter(n, 4, 1);
  }

  /* general case, first find the minimal set of vertices that
   * contain the wiggly lines */
  g0 = dg_mintop(g);
  gc = dg_open(g0->n); /* for fixed edges, empty initially */
  sc = (par ? -1 : 1) + dg_rhsc_recur(g0, gc, par, lookup);
  /* use the Ree-Hoover formula to go from n0 to n */
  sc = dg_iter(n, g0->n, sc);
  dg_close(g0);
  dg_close(gc);
  return sc;
}



/* compute the star content by a look up table */
INLINE int dg_rhsc_lookup(const dg_t *g)
{
  int k, n = g->n;
  dgmap_t *m = dgmap_ + n;
  code_t c;
  static int *sc[DGMAP_NMAX + 1] = {NULL}; /* biconnectivity of unique diagrams */

  if (sc[n] == NULL) { /* initialize the look-up table */
    dg_t *g1;
    int cnt = 0, nz = 0;
    clock_t t0 = clock();

    if (n >= 8) printf("n %d: initializing...\n", n);
    dgmap_init(m, n); /* compute the permutation mapping */
    if (sc[n] == NULL) xnew(sc[n], m->ng);

    if (n >= 8) printf("n %d: diagram-map initialized %gs\n", n, 1.*(clock() - t0)/CLOCKS_PER_SEC);
    /* loop over unique diagrams */
    g1 = dg_open(n);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g1, &m->first[k]);
      if ( dg_biconnected_lookup(g1) ) { /* use the look-up version */
        sc[n][k] = dg_rhsc_low(g1, 1);
        cnt++;
        nz += (sc[n][k] != 0);
      } else sc[n][k] = 0;
    }
    printf("n %d, computed star contents of %d/%d biconnected diagrams\n", n, cnt, nz);
    dg_close(g1);
  }
  dg_encode(g, &c);
  return sc[n][ m->map[c] ];
}



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed */
static int dg_rhsc(dg_t *g)
{
  return (g->n <= DGMAP_NMAX) ? dg_rhsc_lookup(g) : dg_rhsc_low(g, 0);
}



#endif
