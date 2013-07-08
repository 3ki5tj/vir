#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content */



#include "dg.h"
#include "dgdb.h"



/* recursively find the star content */
static int dg_rhsclow(const dg_t *a, int par, int level, dgls_t *ls)
{
  int i, j, sc = 0, add = 0, n = a->n;
  dg_t *b = dg_clone(a);

  for (i = 0; i < n; i++) {
    /* if the degree <= 2, removing the an edge connecting i
     * makes the diagram impossible to biconnected */
    if (dg_deg(b, i) <= 2) continue;
    for (j = i + 1; j < n; j++) {
      if (!dg_linked(b, i, j) || dg_deg(b, j) <= 2)
        continue;
      dg_unlink(b, i, j);
      if ( dg_biconnected(b) ) {
        add = 1;
        if (level >= 1) {
          /* limit the look-up table only to diagrams with the
           * the same number of edges */
          if (dgls_find(ls + level, b->c, b->n) < 0)
            dgls_add(ls + level, b->c, b->n);
          else add = 0;
        }
        if (add) sc += !par ? -1 : 1;
/*
        if (add) {
          printf("i %d, j %d, sc %d, par %d, level %d, ls->n %d\n",
              i, j, sc, par, level, ls[level].n);
          //dg_print(b); getchar();
        }
*/
        sc += dg_rhsclow(b, !par, level + 1, ls);
      }
      dg_link(b, i, j);
    }
  }
  dg_close(b);
  return sc;
}



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



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed */
static int dg_rhsc0(dg_t *g)
{
  int par = 0, sc, i, n = g->n, npr, ned;
  dgls_t *ls;
  dg_t *g0 = NULL;

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

  npr = g0->n * (g0->n - 1)/2;
  xnew(ls, npr); /* list of visited diagrams */
  sc = (par ? -1 : 1) + dg_rhsclow(g0, par, 0, ls);
  for (i = 0; i < npr; i++)
    if (ls->code) free(ls->code);
  free(ls);
  /* use the Ree-Hoover formula to go from n0 to n */
  sc = dg_iter(n, g0->n, sc);
  dg_close(g0);
  return sc;
}



/* compute the star content by a look up table */
INLINE int dg_rhsc_lookup(const dg_t *g)
{
  int k, n = g->n, cnt;
  code_t c;
  dgmap_t *m = dgmap_ + n;
  static int *sc[DGMAP_MAX + 1] = {NULL}; /* biconnectivity of unique diagrams */

  if (sc[n] == NULL) {
    dg_t *g1 = dg_open(n);
    int cnt = 0;

    dgmap_init(m, n); /* compute unique maps */
    xnew(sc[n], m->ng);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g1, &m->first[k]);
      if ( dg_biconnected(g1) ) {
        //dg_print(g1); printf("k %d\n", k);
        sc[n][k] = dg_rhsc0(g1);
        cnt++;
      } else sc[n][k] = 0;
      //dg_print(g1); printf("n %d, k %d, bc %d sc %d\n", n, k, dg_biconnected(g1), sc[n][k]); getchar();
    }
    printf("n %d, computed star contents of %d biconnected diagrams\n", n, cnt);
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
  if (g->n <= DGMAP_MAX) return dg_rhsc_lookup(g);
  else return dg_rhsc0(g);
}



#endif
