#ifndef DGCSEP_H__
#define DGCSEP_H__
/* find clique separator */



#include "dgmap.h"



/* compute a minimal order and the corresponding fill-in of a graph
 * the minimal order is an order of removing vertices such than the
 * resulting edges by contraction is locally a minimum.
 * the algorithm used here constructs a hierarchical spanning tree
 * ``Algorithmic aspects of vertex elimination on graphs''
 * Donald J. Rose, R. Endre Tarjan and George S. Lueker
 * SIAM J. Comput. Vol 5, No. 2, June 1976 */
INLINE void dg_minimalorder(const dg_t *g, dg_t *f, int *a)
{
  int i, j, k, v = 0, w, z, l;
  dgword_t numbered, reached, bv, bw, bz, r;
  DG_DEFN_(g);
  DG_DEFMASKN_();
  int l2[DG_NMAX]; /* the label l times 2 */
  dgword_t reach[DG_NMAX]; /* reach[label] gives a set of vertices */
  int cnt[DG_NMAX * 2];

  if (f) dg_copy(f, g);
  for (i = 0; i < DG_N_; i++) l2[i] = 0; /* l(i) * 2 */
  reached = numbered = 0;
  k = 1; /* number of different labels */
  for (i = DG_N_ - 1; i >= 0; i--) {
    /* select the vertex with the largest label */
    bw = numbered;
    for (r = ~numbered & DG_MASKN_; r; r ^= bv) {
      BITFIRSTLOW(v, r, bv);
      if ( l2[v]/2 == k - 1 ) { /* found it */
        numbered |= bv;
        break;
      }
    }
    die_if (bw == numbered, "i %d no vertex is left\n", i);
    if (a) a[i] = v;
    /* reach[l] gives the set of vertcies with the same label `l'
     * the label `l' is the hierarchy level of the spanning tree
     * vertices with the same label are treated as the same */
    for (l = 0; l < k; l++) reach[l] = 0;
    /* set all unnumbered vertices as unreached */
    reached = numbered;

    /* handle immediate neighbors of v */
    for (r = g->c[v] & ~numbered & DG_MASKN_; r; r ^= bw) {
      BITFIRSTLOW(w, r, bw);
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
        BITFIRSTLOW(w, reach[j], bw);
        reach[j] ^= bw; /* remove w from reach[j] */
        while ( (r = (g->c[w] & ~reached)) != 0 ) {
          BITFIRSTLOW(z, r, bz);
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
    for (l = 0; l < 2*DG_N_; l++) cnt[l] = 0;
    /* accumulate the number of visits of each label */
    for (r = ~numbered & DG_MASKN_; r; r ^= bw) {
      BITFIRSTLOW(w, r, bw);
      cnt[ l2[w] ] = 1;
    }
    /* count the number of different labels */
    for (k = 0, l = 0; l < 2*DG_N_; l++)
      if ( cnt[l] ) /* compute the new label */
        cnt[l] = k++; /* cnt[l] is now the new label */
    for (r = ~numbered & DG_MASKN_; r; r ^= bw) {
      BITFIRSTLOW(w, r, bw);
      l2[w] = 2 * cnt[ l2[w] ]; /* set the new label */
    }
  }
}



/* decompose a diagram by clique separators
 * A clique separator is a fully-connected subgraph
 * `g' is the input diagram, `f' is the fill-in diagram
 * `a' is the elimination order, a[0] is the first vertex to eliminate
 * return the number of cliques, `cl' is the array of cliques
 * `stop1' means stop the search after the first clique separator
 * The algorithm first find a minimal order of elimination.
 * Using this order on a graph with a clique separator, at least
 * one part of the graph is eliminated before the clique
 * ``Decomposition by clique separators'' Robert E. Tarjan,
 * Discrete Mathematics 55 (1985) 221-232 */
INLINE int dg_decompcliqueseplow(const dg_t *g, const dg_t *f,
    const int *a, dgword_t * RESTRICT cl, int stop1)
{
  int v, w, i, ncl = 0;
  DG_DEFN_(g);
  DG_DEFMASKN_();
  dgword_t cb, c, bw, r, unvisited = DG_MASKN_;

  for (i = 0; i < DG_N_; i++) {
    v = a[i];
    unvisited ^= MKBIT(v); /* remove the `v' bit */
    /* compute C(v), the set of succeeding vertices that
     * are adjacent to v */
    c = unvisited & f->c[v];
    /* test if C(v) is a clique, a fully-connected subgraph */
    for (r = c; r; r ^= bw) {
      BITFIRSTLOW(w, r, bw);
      /* c ^ bw is the set of vertices connected to `w'
       * in `c', if `c' is a clique */
      cb = c ^ bw;
      if ((g->c[w] & cb) != cb) /* not a clique */
        break; /* break the loop prematurally, r != 0 */
    }
    if (r == 0) { /* if the loop is completed, `c' is a clique */
      if (unvisited  == c) { /* clique `c' == the rest vertices */
        return ncl;          /* so it is not a separator */
      } else { /* found a clique `c' */
        if (cl != NULL) cl[ncl] = c;
        ncl++;
        if (stop1) return 1;
      }
    }
  }
  return ncl;
}



/* test if a graph has a clique separator */
INLINE dgword_t dg_cliquesep(const dg_t *g)
{
  static dg_t *fs[DG_NMAX + 1];
  static int a[DG_NMAX]; /* a[k] is the kth vertex */
#pragma omp threadprivate(fs, a)
  dg_t *f; /* fill-in graph */
  DG_DEFN_(g);
  dgword_t cl;

  if (fs[DG_N_] == NULL) fs[DG_N_] = dg_open(DG_N_);
  f = fs[DG_N_];

  /* 1. find a minimal ordering and its fill-in */
  dg_minimalorder(g, f, a);

  /* 2. clique decomposition (stop after the first clique) */
  if ( dg_decompcliqueseplow(g, f, a, &cl, 1) ) return cl;
  else return 0;
}



/* number of nodes in the clique-separator decomposition */
#define dg_ncsep(g) dg_decompcsep(g, NULL)

INLINE int dg_decompcsep(const dg_t *g, dgword_t * RESTRICT cl)
{
  static dg_t *fs[DG_NMAX + 1];
  static int a[DG_NMAX]; /* a[k] is the kth vertex */
#pragma omp threadprivate(fs, a)
  dg_t *f;
  DG_DEFN_(g);

  if (fs[DG_N_] == NULL) fs[DG_N_] = dg_open(DG_N_);
  f = fs[DG_N_];

  /* 1. find a minimal ordering and its fill-in */
  dg_minimalorder(g, f, a);

  /* 2. clique decomposition */
  return dg_decompcliqueseplow(g, f, a, cl, 0);
}



#ifdef DGMAP_EXISTS
/* compute the number of nodes the clique-separator decomposition */
INLINE int dg_ncsep_lookuplow(const dg_t *g, dgword_t c)
{
  static char *ncl[DGMAP_NMAX + 1];
#pragma omp threadprivate(ncl)
  int ncs;
  DG_DEFN_(g);

  /* initialize the lookup table */
  if (ncl[DG_N_] == NULL) {
    int ipr, npr = 1u << (DG_N_ * (DG_N_ - 1) / 2);
    xnew(ncl[DG_N_], npr);
    for (ipr = 0; ipr < npr; ipr++) ncl[DG_N_][ipr] = (char) (-1);
  }

  ncs = ncl[DG_N_][c];
  if (ncs < 0) {
    ncs = dg_ncsep(g);
    ncl[DG_N_][c] = (char) ncs; /* save the value */
  }
  return ncs;
}


/* compute the number of nodes the clique-separator decomposition */
INLINE int dg_ncsep_lookup(const dg_t *g)
{
  dgword_t code;

  die_if (g->n > DGMAP_NMAX, "n %d too large\n", g->n);
  dg_encode(g, &code);
  return dg_ncsep_lookuplow(g, code);
}
#endif /* defined(DGMAP_EXISTS) */

#endif

