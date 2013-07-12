#ifndef DGCSEP_H__
#define DGCSEP_H__
/* find clique separator */

#include "dg.h"


/* compute a minimal order */
INLINE void dg_minimalorder(const dg_t *g, dg_t *f, int *a, int *p)
{
  int i, j, k, v, w, z, l, n = g->n;
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
      if (f) dg_link(f, v, w);
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
            if (f) dg_link(f, v, z);
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



/* test if a graph has a clique separator */
INLINE code_t dg_cliquesep(const dg_t *g)
{
  static dg_t *fs[DG_NMAX + 1], *f;
  static int a[DG_NMAX], p[DG_NMAX]; /* a[k] is the kth vertex, p = a^(-1) */
  int i, v, w, n = g->n;
  code_t c, bw, r, unvisited;

  if (fs[n] == NULL) f = fs[n] = dg_open(n);

  /* 1. find a minimal ordering and its fill-in */
  dg_minimalorder(g, f, a, p);

  /* 2. clique decomposition (stop after the first clique) */
  for (unvisited = ((code_t) 1 << n) - 1, i = 0; i < n; i++) {
    v = a[i];
    unvisited &= ~((code_t) 1 << v);
    /* compute C(v) */
    for (c = 0, w = 0; w < n; w++)
      if (dg_linked(f, v, w) && p[w] > p[v])
        c |= 1 << w;
    /* test if C(v) is a clique, a fully-connected subgraph */
    for (r = c; r; r ^= bw) {
      w = bitfirstlow(r, &bw);
      if ((g->c[w] & (bw ^ c)) != (bw ^ c))
        break;
    }
    if (r == 0) {
      if (unvisited == c) return 0; /* the rest */
      else return c; /* found a clique */
    }
  }
  return 0;
}



#endif

