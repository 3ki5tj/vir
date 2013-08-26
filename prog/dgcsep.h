#ifndef DGCSEP_H__
#define DGCSEP_H__
/* find clique separator */

#include "dg.h"


/* compute a minimal order and the corresponding fill-in of a graph
 * the minimal order is an order of removing vertices such than the
 * resulting edges by contraction is locally a minimum.
 * the algorithm used here constructs a hierarchical spanning tree
 * ``Algorithmic aspects of vertex elimination on graphs''
 * Donald J. Rose, R. Endre Tarjan and George S. Lueker
 * SIAM J. Comput. Vol 5, No. 2, June 1976 */
INLINE void dg_minimalorder(const dg_t *g, dg_t *f, int *a, int *p)
{
  int i, j, k, v = 0, w, z, l, n = g->n;
  code_t numbered, reached, bv, bw, bz, r, mask;
  int l2[DG_NMAX]; /* the label l times 2 */
  code_t reach[DG_NMAX]; /* reach[label] gives a set of vertices */
  int cnt[DG_NMAX * 2];

  mask = ((code_t) 1 << n) - 1;
  if (f) dg_copy(f, g);
  for (i = 0; i < n; i++) l2[i] = 0; /* l(i) * 2 */
  reached = numbered = 0;
  k = 1; /* number of different labels */
  for (i = n - 1; i >= 0; i--) {
    /* select the vertex with the largest label */
    bw = numbered;
    for (r = ~numbered & mask; r; r ^= bv) {
      v = bitfirstlow(r, &bv);
      if ( l2[v]/2 == k - 1 ) { /* found it */
        numbered |= bv;
        break;
      }
    }
    die_if (bw == numbered, "i %d no vertex is left\n", i);
    if (a) a[i] = v;
    if (p) p[v] = i;
    /* reach[l] gives the set of vertcies with the same label `l'
     * the label `l' is the hierarchy level of the spanning tree
     * vertices with the same label are treated as the same */
    for (l = 0; l < k; l++) reach[l] = 0;
    /* set all unnumbered vertices as unreached */
    reached = numbered;

    /* handle immediate neighbors of v */
    for (r = g->c[v] & ~numbered & mask; r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      /* the immediate neighbors of w are the starting points
       * of the search, we group them by the labels */
      reach[ l2[w]/2 ] |= bw;
      reached |= bw;
      /* only update the label after we have done with it
       * the l2[w]++ operation only applies to a vertex once */
      l2[w]++;
      if (f) DG_LINK(f, v, w);
    }

    /* search paths from v
     * the loop is over all labels, from smaller labels to larger ones.
     * during the search, only vertices of equal or larger labels are
     * produced, so the one-pass search is sufficient */
    for (j = 0; j < k; j++) {
      while ( reach[j] ) { /* if the vertex set is not empty */
        /* delete the first w from reach[j] */
        w = bitfirstlow(reach[j], &bw);
        reach[j] ^= bw; /* remove w from reach[j] */
        while ( (r = (g->c[w] & ~reached)) != 0 ) {
          z = bitfirstlow(r, &bz);
          reached |= bz;
          if ((l = l2[z]/2) > j) {
            reach[l] |= bz;
            l2[z]++;
            if (f) DG_LINK(f, v, z);
          } else { /* lower label encountered, count it as label j */
            reach[j] |= bz;
          }
        }
      }
    }

    /* re-assign labels of the vertices by counting sort */
    for (l = 0; l < 2*n; l++) cnt[l] = 0;
    /* accumulate the number of visits of each label */
    for (r = ~numbered & mask; r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      cnt[ l2[w] ] = 1;
    }
    /* count the number of different labels */
    for (k = 0, l = 0; l < 2*n; l++)
      if ( cnt[l] ) /* compute the new label */
        cnt[l] = k++; /* cnt[l] is now the new label */
    for (r = ~numbered & mask; r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      l2[w] = 2 * cnt[ l2[w] ]; /* set the new label */
    }
  }
}



/* decompose a diagram by clique separators
 * A clique separator is a fully-connected subgraph
 * `g' is the input diagram, `f' is the fill-in diagram
 * return the number of cliques, `cl' is the array of cliques
 * `stop1' means stop the search after the first step
 * The algorithm first find a minimal order of elimination.
 * Using this order on a graph with a clique separator, at least
 * one part of the graph is eliminated before the clique
 * ``Decomposition by clique separators'' Robert E. Tarjan,
 * Discrete Mathematics 55 (1985) 221-232 */
INLINE int dg_decompcliqueseplow(const dg_t *g, const dg_t *f,
    const int *a, const int *p, code_t * RESTRICT cl, int stop1)
{
  int v, w, i, n = g->n, ncl = 0;
  code_t cb, c, bw, r, unvisited;

  for (unvisited = ((code_t) 1u << n) - 1, i = 0; i < n; i++) {
    v = a[i];
    unvisited &= ~((code_t) 1 << v); /* remove the `v' bit */
    /* compute C(v), the set of succeeding vertices that
     * are adjacent to v */
    for (c = 0, r = (unvisited & f->c[v]); r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      if (p[w] > p[v])
        c |= (code_t) 1u << w;
    }
    /* test if C(v) is a clique, a fully-connected subgraph */
    for (r = c; r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      /* c ^ bw is the set of vertices connected to `w'
       * in `c', if `c' is a clique */
      cb = c ^ bw;
      if ((g->c[w] & cb) != cb) /* not a clique */
        break; /* break the loop prematurally */
    }
    if (r == 0) { /* if the loop completed, `c' is a clique */
      if (unvisited  == c) { /* clique `c' == the rest vertices */
        return ncl;
      } else { /* found a clique `c' */
        if (cl) cl[ncl] = c;
        ncl++;
        if (stop1) return 1;
      }
    }
  }
  return ncl;
}



/* test if a graph has a clique separator */
INLINE code_t dg_cliquesep(const dg_t *g)
{
  static dg_t *fs[DG_NMAX + 1], *f; /* fill-in graph */
  static int a[DG_NMAX], p[DG_NMAX]; /* a[k] is the kth vertex, p = a^(-1) */
  int n = g->n;
  code_t cl;

  if (fs[n] == NULL) fs[n] = dg_open(n);
  f = fs[n];

  /* 1. find a minimal ordering and its fill-in */
  dg_minimalorder(g, f, a, p);

  /* 2. clique decomposition (stop after the first clique) */
  if ( dg_decompcliqueseplow(g, f, a, p, &cl, 1) ) return cl;
  else return 0;
}



/* number of nodes in the clique-separator decomposition */
#define dg_ncsep(g) dg_decompcsep(g, NULL)

INLINE int dg_decompcsep(const dg_t *g, code_t * RESTRICT cl)
{
  static dg_t *fs[DG_NMAX + 1], *f;
  static int a[DG_NMAX], p[DG_NMAX]; /* a[k] is the kth vertex, p = a^(-1) */
  int n = g->n;

  if (fs[n] == NULL) fs[n] = dg_open(n);
  f = fs[n];

  /* 1. find a minimal ordering and its fill-in */
  dg_minimalorder(g, f, a, p);

  /* 2. clique decomposition */
  return dg_decompcliqueseplow(g, f, a, p, cl, 0);
}



/* compute the number of nodes the clique-separator decomposition */
INLINE int dg_ncsep_lookuplow(const dg_t *g, code_t c)
{
  static char *ncl[DGMAP_NMAX + 1];
  int n = g->n;

  if (ncl[n] == NULL) {
    int ipr, npr = (code_t) 1u << (n * (n - 1) / 2);
    xnew(ncl[n], npr);
    for (ipr = 0; ipr < npr; ipr++) ncl[n][ipr] = (char) (-1);
  }
  if (ncl[n][c] < 0) ncl[n][c] = dg_ncsep(g);
  return ncl[n][c];
}



/* compute the number of nodes the clique-separator decomposition */
INLINE int dg_ncsep_lookup(const dg_t *g)
{
  code_t code;

  die_if (g->n > DGMAP_NMAX, "n %d too large\n", g->n);
  dg_encode(g, &code);
  return dg_ncsep_lookuplow(g, code);
}

#endif

