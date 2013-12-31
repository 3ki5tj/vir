#ifndef DGCSEP_H__
#define DGCSEP_H__



/* find clique separators */



#include "dgmap.h"



/* methods of testing clique separators */
enum {
  DGCSEP_NULLMETHOD = 0, /* disable the detection of clique separators */
  DGCSEP_LEXM = 1, /* lexicographic search for minimal ordering */
  DGCSEP_MCSP = 2, /* maximal cardinality search for perfect ordering */
  DGCSEP_MCSM = 3, /* maximal cardinality search for minimal ordering */
  DGCSEP_NBASICMETHODS = 4
};

#ifndef DGCSEP_DEFAULTMETHOD
#define DGCSEP_DEFAULTMETHOD 1
#endif

/* ways of finding clique separators during the computation of minimal ordering */
enum {
  DGCSEP_NOFIND = 0,
  DGCSEP_FIND1 = 1,
  DGCSEP_FINDALL = 2
};



/* check if vertex v is adjacent to a clique among the set `vs'
 * this function embeds the second step of Tarjan's algorithm
 * into a routine of finding minimal ordering
 * ``Decomposition by clique separators'' Robert E. Tarjan,
 * Discrete Mathematics 55 (1985) 221-232 */
INLINE dgvsref_t dg_adjclique(const dg_t *g, const dg_t *f, int v, dgword_t vs)
{
  static dgvs_t c;
#pragma omp threadprivate(c)
  dgvs_t r, cb, cb1;
  dgword_t bz;
  int z;
  DGVS_DEFIQ_(zq)

  DGVS_AND2(c, f->c[v], vs)

  /* clique `c' == the entire vertex set,
   * `c' is not a separator even if it is a clique */
  if ( dgvs_eq(c, vs) ) return 0;

  /*  test if c is a clique separator */
  DGVS_CPY(r, c)
  DGVS_CPY(cb, c)
  while ( dgvs_nonzero(r) ) {
    DGVS_FIRSTLOW(z, r, bz, zq)
    DGVS_XOR1(r, bz, zq) /* remove vertex `z' */
    /* c ^ bz is the set of vertices connected to `z' in `c',
     * if `c' is a clique */
    DGVS_XOR1(cb, bz, zq)
    DGVS_AND2(cb1, g->c[z], cb)
    if ( dgvs_neq(cb1, cb) ) /* not a clique */
      return 0; /* break the loop prematurely, r != 0 */
    DGVS_XOR1(cb, bz, zq)
  }
  /* if the loop is completed, `c' is a clique */
  return c;
}



/* maximal cardinality search
 * Robert E. Tarjan and Mihalis Yannakakis
 * ``Simple Linear-Time Algorithms to Test Chordality of
 *   Graphs, Test Acyclicity of Hypergraphs, and Selectively
 *   Reduce Acyclic Hypergraphs''
 * SIAM J. Comput. Vol. 13 No. 3 August 1984 */
INLINE void dg_mcsearchp(const dg_t *g, int *a, int *p)
{
  static dgvs_t set[DG_NMAX];
  static int size[DG_NMAX];
#pragma omp threadprivate(set, size)
  int i, j, v, w, szw;
  dgvs_t unvisited;
  dgword_t cv, bv, bw;
  DGVS_DEFIQ_(vq)
  DGVS_DEFIQ_(wq)
  DG_DEFN_(g)

  /* `set[i]' collects vertices of the same cardinality,
   * which is the number of neighboring numbered vertices */
  for (i = 0; i < DG_N_; i++) {
    DGVS_CLEAR( set[i] )
  }
  /* `size[v]' gives the cardinality of vertex `v' */
  for (i = 0; i < DG_N_; i++) size[i] = 0;
  DGVS_MKBITSMASK(set[0], DG_N_)
  DGVS_MKBITSMASK(unvisited, DG_N_)
  j = 0;
  for (i = DG_N_ - 1; i >= 0; i--) { /* number the ith vertex */
    /* set[j] collects vertices with the highest cardinality */
    DGVS_FIRSTLOW(v, set[j], bv, vq)
    DGVS_XOR1(set[j], bv, vq)
    DGVS_XOR1(unvisited, bv, vq)
    a[i] = v;
    p[v] = i;
    DGVS_AND2(cv, g->c[v], unvisited)
    while ( dgvs_nonzero(cv) ) {
      DGVS_FIRSTLOW(w, cv, bw, wq)
      DGVS_XOR1(cv, bw, wq)
      szw = size[w];
      DGVS_XOR1(set[szw], bw, wq) /* delete `w' from `set[szw]' */
      DGVS_XOR1(set[size[w] = szw + 1], bw, wq) /* add `w' to `set[szw+1]' */
    }
    /* seek the next largest cardinality */
    for (j++; j >= 0 && set[j] == 0; j--)
      ;
  }
}



/* compute the fill-in graph `f' according to the ordering `a' */
INLINE void dg_fillin(const dg_t *g, dg_t *f,
    const int *a, const int *p)
{
  static int index[DG_NMAX], follow[DG_NMAX];
#pragma omp threadprivate(index, follow)
  int i, w, v, x;
  dgword_t cw, bv;
  DGVS_DEFIQ_(vq)
  DG_DEFN_(g)

  dg_empty(f);
  for (i = 0; i < DG_N_; i++) {
    w = a[i];
    /* `follow[w]' gives the first succeeding neighbor of `w'
     * in the fill-in graph `f' */
    follow[w] = w;
    /* `index[v]' gives the order of the last adjacent neighbor
     * of `v' in the original graph `g' */
    index[w] = i;
    DGVS_CPY(cw, g->c[w])
    while ( dgvs_nonzero(cw) ) {
      DGVS_FIRSTLOW(v, cw, bv, vq)
      DGVS_XOR1(cw, bv, vq)
      if (p[v] >= i) continue;
      for (x = v; index[x] < i; x = follow[x] ) {
        index[x] = i;
        DG_LINK(f, x, w)
      }
      if (follow[x] == x)
        follow[x] = w;
    }
  }
}



/* compute a minimal ordering and the corresponding elimination graph
 * by algorithm MCS-M (maximal cardinality search)
 * the minimal ordering is an ordering of removing vertices such than the
 * resulting edges by contraction is marginally minimized
 * Anne Berry, Jean R. S. Blair and Pinar Heggernes
 * ``Maximum Cardinality Search for Computing Minimal Triangulations''
 *   in Graph-Theoretic Concepts in Computer Science,
 *   Vol. 2573 pg. 1-12 (2002).
 * If `csepmode' is 1, test clique separators and stop at the first find
 *               is 2, find all clique separators in `cl'
 * return the number of clique separators found */
INLINE int dg_minimalordermcsm(const dg_t *g, dg_t *f, int *a,
    int csepmode, dgword_t *cl)
{
  static dgvs_t set[DG_NMAX];
  static int size[DG_NMAX];
  static dgvs_t reach[DG_NMAX]; /* reach[size] gives a queue of vertices
                                     of the same cutoff cardinality size */
#pragma omp threadprivate(set, size, reach)
  int i, j, k, v = 0, w, z, l, szw;
  dgvs_t unnumbered, reached, r;
  dgvsref_t reachj;
  dgword_t bv, bw, bz;
  int ncl = 0; /* number of clique separators */
  DGVS_DEFIQ_(vq)
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(iq1)
  DG_DEFN_(g)
  DGVS_DEFMASKN_()

  dg_copy(f, g);
  /* `set[i]' collects vertices of the same cardinality,
   * which is the number of neighboring numbered vertices */
  for (i = 0; i < DG_N_; i++) {
    DGVS_CLEAR( set[i] )
  }
  DGVS_CPY(set[0], DGVS_MASKN_)
  /* `size[v]' gives the cardinality of vertex `v' */
  for (i = 0; i < DG_N_; i++) size[i] = 0;
  k = 0; /* maximal cardinality at the moment */
  i = DG_N_ - 1;
  DGVS_CPY(unnumbered, DGVS_MASKN_)
  while (1) {
#ifdef DGCSEP_DBG
    dgvs_t unnumbered0;
    DGVS_CPY(unnumbered0, unnumbered)
#endif
    reachj = set[k];
    DGVS_FIRSTLOW(v, reachj, bv, vq)
    DGVS_XOR1(reachj, bv, vq)
    DGVS_CPY(set[k], reachj)
    a[i] = v;

    /* test clique separators */
    if ( csepmode != 0 ) {
      dgvs_t numbered;
      dgvsref_t clq;
      DGVS_XOR2(numbered, unnumbered, DGVS_MASKN_)
      if ( (clq = dg_adjclique(g, f, v, numbered)) != 0 ) {
        if ( cl != NULL ) {
          DGVS_CPY(cl[ncl], clq)
        }
        if ( csepmode == 1 ) return 1;
        ncl++;
      }
    } /* (cl != NULL) */

    DGVS_REMOVE(unnumbered, v)

    /* I. Preparation */
    /* reach[l] gives the set of vertices with the same size `l'
     * the size `l' is the hierarchy level of the spanning tree
     * Vertices with the same size are treated as the equivalent,
     * and thus are grouped together
     * In block III below, we loop over levels instead of vertices  */
    for (l = 0; l < k; l++)
      DGVS_CLEAR( reach[l] )
    /* set all unnumbered vertices as unreached */
    DGVS_XOR2(reached, unnumbered, DGVS_MASKN_)

    /* II. Handle unnumbered vertices `v' adjacent to `v' in the graph `g'
     * The above unnumbered neighbors are the seeds in constructing
     * the hierarchical spanning tree based on the breath-first search */
    DGVS_AND2(r, g->c[v], unnumbered) /* r = g->c[v] & unnumbered */
    /* mark vertices in `r' as `reached' such that the size of
     * these wil not be updated further in this round */
    DGVS_OR(reached, r) /* reached |= r */
    while ( dgvs_nonzero( r ) ) { /* loop over vertices `w' in `r' */
      DGVS_FIRSTLOW(w, r, bw, iq) /* first unnumbered vertex `w' */
      DGVS_XOR1(r, bw, iq) /* remove the vertex `w' from the set `r' */
      szw = size[w];
      /* the immediate neighbors of `w' are the starting points
       * of the search below, we group them by the size */
      DGVS_OR1(reach[ szw ], bw, iq)
      /* increase the size of `w' such that vertices reached by `v'
       * will have larger sizes than others
       * At the end of this loop, there might be multiple vertices
       * sharing the same odd size; but this is okay, for these
       * vertices of the same size are equivalent, i.e., they
       * are all adjacent to `v' in the fill-in graph */
      DGVS_XOR1(set[szw], bw, iq)
      size[w] = ++szw;
      DGVS_XOR1(set[szw], bw, iq)
      /* add the edge (`v', `w') into the fill-in graph `f' */
      DG_LINK(f, v, w);
    }

    /* III. Search unnumbered indirect neighbors of `v', which are vertices
     * adjacent to `v' in the fill-in graph `f', but not in the original
     * graph `g'
     * We will search paths v - z1 (= w) - z2 - ... - zk - z{k+1} (= z),
     * such that l2(z1), ..., l2(zk) are all greater than l2(z{k+1}),
     * and z1 is one of * the seeds `w' found in the previous block (block II)
     * cf Lemma 7 and Sec. 4.1 on page 273 of the RTL paper
     * The search over such paths is implemented as a BFS rooted at
     * z1 = w, which is called the seed vertex.
     * This allows us to reach all intermediate vertices with sizes
     * no larger than the size of w, which is j, until the terminal
     * vertex v is found.
     * To avoid repeated searches over the same path segments, we need to
     * first sort the seed vertex z1 in ascending order of the size.
     * This sets the lowest threshold for subsequent vertices on the paths.
     * Once a vertex is reached from the seed vertex with the smallest
     * size, it is unnecessary to search paths from vertices of
     * larger sizes.
     * Instead of being over all vertices, the loop is over all sizes,
     * from smaller sizes to larger ones.  This is essentially a radix
     * sort of the seed vertices.
     * The updates `l2[w]++' in the previous block (block II) are irrelevant
     * here, because they happen after the array reach[] is updated
     * */
    for (j = 0; j < k; j++) { /* loop over distinct labels */
      /* loop over vertices with the same label j
       * note, the vertex set reach[j] will be updated within the loop
       * so it should not be fixed before the loop */
      reachj = reach[j]; /* make a reference point */
      while ( dgvs_nonzero( reachj ) ) {
        /* pop the first `w' from reachj */
        DGVS_FIRSTLOW(w, reachj, bw, iq)
        DGVS_XOR1(reachj, bw, iq) /* remove w from reachj */
        /* loop over unreached vertices adjacent to g->c[w] */
        DGVS_ANDNOT2(r, g->c[w], reached) /* r = g->c[w] & ~reached */
        DGVS_OR(reached, r) /* mark vertices adjacent to `w' as reached */
        while ( dgvs_nonzero( r ) ) {
          DGVS_FIRSTLOW(z, r, bz, iq1) /* pick `z' */
          DGVS_XOR1(r, bz, iq1) /* remove `z' from `r' */
          /* here we distinguish two cases
           * (A) if `z' has a label greater than `j' (the label of `w')
           * it means that the `w' is eliminated before `z', so there
           * will be an edge `v'-`z' in the fill-in graph
           * Thus, we should add vertex `z' (just as we have added `w')
           * to the array reach[] for further searches
           * Since `z' has a larger label, its label will be handled later
           * in the outer loop
           * (B) if, however, `z' has a label no larger than `j' (the label of `w')
           * it means that `z' is eliminated before `w'
           * But instead of going back to the smaller label of `z',
           * we simply append `z' to this list reach[j] */
          if ((l = size[z]) > j) {
            DGVS_OR1(reach[l], bz, iq1)
            DGVS_XOR1(set[l], bz, iq1)
            size[z] = l + 1;
            DGVS_XOR1(set[l + 1], bz, iq1)
            DG_LINK(f, v, z);
          } else { /* lower label encountered, count it as label[j] */
            DGVS_OR1(reachj, bz, iq1)
          }
        }
      }
    }

    if (i == 0) break;
    i--;

    /* find the set with the largest cardinality */
    for (k++; k >= 0 && dgvs_iszero(set[k]); k--)
      ;
  }

  return ncl;
}



/* rearrange l2[] such that they occupied 0..k-1 positions */
INLINE int dg_sortlabels(int l2[], dgvs_t unnumbered, int n)
{
  int l, w, k, lmax = 0;
  int cnt[DG_NMAX * 2];
  dgvs_t r;
  dgword_t bw;
  DGVS_DEFIQ_(iq)

  /* A. Check if each label is visited */
  for (l = 0; l < n*2; l++) cnt[l] = 0;
  DGVS_CPY(r, unnumbered)
  while ( dgvs_nonzero(r) ) { /* loop over unnumbered vertices in `r' */
    DGVS_FIRSTLOW(w, r, bw, iq) /* extract the first vertex `w' in `r' */
    l = l2[w];
    if (l > lmax) lmax = l;
    cnt[ l ] = 1;
    DGVS_XOR1(r, bw, iq) /* remove `w' from the vertex set `r' */
  }
  /* B. Count the number of different labels */
  for (k = 0, l = 0; l <= lmax; l++)
    if ( cnt[l] ) /* compute the new label */
      cnt[l] = k++; /* cnt[l] is now the new label */
#ifdef DGCSEP_DBG
  printf("counting sort %d items, n %d, number of values %d, lmax %d\n",
     dgvs_count(unnumbered), n, k, lmax);
#endif
  /* C. Assign the new labels to the array l2[] */
  DGVS_CPY(r, unnumbered)
  while ( dgvs_nonzero(r) ) { /* loop over unnumbered vertices in `r' */
    DGVS_FIRSTLOW(w, r, bw, iq) /* extract the first vertex `w' in `r' */
    l2[w] = 2 * cnt[ l2[w] ]; /* set the new label */
    DGVS_XOR1(r, bw, iq) /* remove `w' from the vertex set `r' */
  }
  return k;
}



#define dg_minimalorder(g, f, a) dg_minimalorderlexm(g, f, a, 0, NULL)

/* compute a minimal ordering and the corresponding elimination graph
 * by algorithm LEX-M (lexicographic search)
 * the minimal ordering is an ordering of removing vertices such than the
 * resulting edges by contraction is marginally minimized.
 * ``Algorithmic aspects of vertex elimination on graphs''
 * Donald J. Rose, R. Endre Tarjan and George S. Lueker
 * SIAM J. Comput. Vol 5, No. 2, June 1976
 * If `csepmode' is 1, test clique separators and stop at the first find
 *               is 2, find all clique separators in `cl'
 * return the number of clique separators found */
INLINE int dg_minimalorderlexm(const dg_t *g, dg_t *f, int *a,
    int csepmode, dgword_t *cl)
{
  int i, j, k, v = 0, w, z, l;
  dgvs_t unnumbered, reached, r;
  dgvsref_t reachj;
  dgword_t bw, bz;
  int l2[DG_NMAX]; /* the label l times 2 */
  dgvs_t reach[DG_NMAX]; /* reach[label] gives a set of vertices */
  int ncl = 0; /* number of clique separators */
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(iq1)
  DG_DEFN_(g)
  DGVS_DEFMASKN_()

  dg_copy(f, g);
  for (i = 0; i < DG_N_; i++) l2[i] = 0; /* l(i) * 2 */
  k = 1; /* number of different labels */
  i = DG_N_ - 1;
  w = 0;
  DGVS_MKBIT(w, bw, iq)
  DGVS_CPY(unnumbered, DGVS_MASKN_)
  while (1) {
#ifdef DGCSEP_DBG
    dgvs_t unnumbered0;
    DGVS_CPY(unnumbered0, unnumbered)
#endif
    a[i] = v = w;

    /* test clique separators */
    if ( csepmode != 0 ) {
      dgvs_t numbered;
      dgvsref_t clq;
      DGVS_XOR2(numbered, unnumbered, DGVS_MASKN_)
      if ( (clq = dg_adjclique(g, f, v, numbered)) != 0 ) {
        if ( cl != NULL ) {
          DGVS_CPY(cl[ncl], clq)
        }
        if ( csepmode == 1 ) return 1;
        ncl++;
      }
    } /* (cl != NULL) */

    DGVS_REMOVE(unnumbered, v)

    /* I. Preparation */
    /* reach[l] gives the set of vertices with the same label `l'
     * the label `l' is the hierarchy level of the spanning tree
     * Vertices with the same label are treated as the equivalent,
     * and thus are grouped together
     * In block III below, we loop over levels instead of vertices  */
    for (l = 0; l < k; l++)
      DGVS_CLEAR( reach[l] )
    /* set all unnumbered vertices as unreached */
    DGVS_XOR2(reached, unnumbered, DGVS_MASKN_)

    /* II. Handle unnumbered vertices `v' adjacent to `v' in the graph `g'
     * The above unnumbered neighbors are the seeds in constructing
     * the hierarchical spanning tree based on the breath-first search */
    DGVS_ANDNOT2(r, g->c[v], reached) /* r = g->c[v] & ~reached */
    /* mark vertices in `r' as `reached' such that the label of
     * these wil not be updated further in this round */
    DGVS_OR(reached, r) /* reached |= r */
    while ( dgvs_nonzero( r ) ) { /* loop over vertices `w' in `r' */
      DGVS_FIRSTLOW(w, r, bw, iq) /* first unnumbered vertex `w' */
      /* the immediate neighbors of `w' are the starting points
       * of the search below, we group them by the labels
       * l2[w]/2 gives the group id, l2[w] is even at the moment */
      DGVS_OR1(reach[ l2[w] / 2 ], bw, iq)
      /* increase the label of `w' such that vertices reached by `v'
       * will have larger labels than others
       * Since each label is updated only once, the label l2[w] encountered
       * here is always an even number before the update below
       * At the end of this loop, there might be multiple vertices
       * sharing the same odd label; but this is okay, for these
       * vertices of the same odd label are equivalent, i.e., they
       * are all adjacent to `v' in the fill-in graph */
      l2[w]++;
      /* add the edge (`v', `w') into the fill-in graph `f' */
      DG_LINK(f, v, w);
      DGVS_XOR1(r, bw, iq) /* remove the vertex `w' from the set `r' */
    }

    /* III. Search unnumbered indirect neighbors of `v', which are vertices
     * adjacent to `v' in the fill-in graph `f', but not in the original
     * graph `g'
     * We will search paths v - z1 (= w) - z2 - ... - zk - z{k+1} (= z),
     * such that l2(z1), ..., l2(zk) are all greater than l2(z{k+1}),
     * and z1 is one of * the seeds `w' found in the previous block (block II)
     * cf Lemma 7 and Sec. 4.1 on page 273 of the RTL paper
     * The search over such paths is implemented as a BFS rooted at
     * z1 = w, which is called the seed vertex.
     * This allows us to reach all intermediate vertices with labels
     * no larger than the label of w, which is j, until the terminal
     * vertex v is found.
     * To avoid repeated searches over the same path segments, we need to
     * first sort the seed vertex z1 in ascending order of the label.
     * This sets the lowest threshold for subsequent vertices on the paths.
     * Once a vertex is reached from the seed vertex with the smallest
     * label, it is unnecessary to search paths from vertices of
     * larger labels.
     * Instead of being over all vertices, the loop is over all labels,
     * from smaller labels to larger ones.  This is essentially a radix
     * sort of the seed vertices.
     * The updates `l2[w]++' in the previous block (block II) are irrelevant
     * here, because they happen after the array reach[] is updated
     * */
    for (j = 0; j < k; j++) { /* loop over distinct labels */
      /* loop over vertices with the same label j
       * note, the vertex set reach[j] will be updated within the loop
       * so it should not be fixed before the loop */
      reachj = reach[j]; /* make a reference point */
      while ( dgvs_nonzero( reachj ) ) {
        /* pop the first `w' from reachj */
        DGVS_FIRSTLOW(w, reachj, bw, iq)
        DGVS_XOR1(reachj, bw, iq) /* remove w from reachj */
        /* loop over unreached vertices adjacent to g->c[w] */
        DGVS_ANDNOT2(r, g->c[w], reached) /* r = g->c[w] & ~reached */
        DGVS_OR(reached, r) /* mark vertices adjacent to `w' as reached */
        while ( dgvs_nonzero( r ) ) {
          DGVS_FIRSTLOW(z, r, bz, iq1) /* pick `z' */
          DGVS_XOR1(r, bz, iq1) /* remove `z' from `r' */
          /* here we distinguish two cases
           * (A) if `z' has a label greater than `j' (the label of `w')
           * it means that the `w' is eliminated before `z', so there
           * will be an edge `v'-`z' in the fill-in graph
           * Thus, we should add vertex `z' (just as we have added `w')
           * to the array reach[] for further searches
           * Since `z' has a larger label, its label will be handled later
           * in the outer loop
           * (B) if, however, `z' has a label no larger than `j' (the label of `w')
           * it means that `z' is eliminated before `w'
           * But instead of going back to the smaller label of `z',
           * we simply append `z' to this list reach[j] */
          if ((l = l2[z]/2) > j) {
            DGVS_OR1(reach[l], bz, iq1)
            l2[z]++;
            DG_LINK(f, v, z);
          } else { /* lower label encountered, count it as label[j] */
            DGVS_OR1(reachj, bz, iq1)
          }
        }
      }
    }

    if (i == 0) break;
    i--;

    /* IV. Re-assign labels of the vertices by counting sort */
    k = dg_sortlabels(l2, unnumbered, DG_N_);
#ifdef DGCSEP_DBG
    printf("round %d, %d labels\n", i, k);
#endif

    /* V. label the next vertex a(i) */
    DGVS_CPY(r, unnumbered)
    while ( dgvs_nonzero(r) ) {
      DGVS_FIRSTLOW(w, r, bw, iq)
      DGVS_XOR1(r, bw, iq)
      if ( l2[w] == k*2 - 2 )
        break;
    }
  }

  return ncl;
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
    const int *a, dgvs_t * RESTRICT cl, int stop1)
{
  int v, w, i, ncl = 0;
  dgvs_t cb, c, r, cwb, unvisited;
  dgword_t bw;
  DGVS_DEFIQ_(iq)
  DG_DEFN_(g)
  DG_DEFMASKN_()

  DGVS_CPY(unvisited, DG_MASKN_)
  for (i = 0; i < DG_N_; i++) {
    v = a[i];
    DGVS_FLIP(unvisited, v) /* remove the `v' bit */
    /* compute C(v), the set of succeeding vertices that
     * are adjacent to v */
    DGVS_AND2(c, unvisited, f->c[v]) /* c = unvisited & f->c[v] */
    /* test if c = C(v) is a clique, a fully-connected subgraph */
    DGVS_CPY(r, c)
    while ( dgvs_nonzero(r) ) { /* loop over vertices in `c' */
      DGVS_FIRSTLOW(w, r, bw, iq) /* vertex `w' in `c' */
      /* cb = c ^ bw is the set of vertices connected to `w' in `c'
       * if `c' is a clique */
      DGVS_CPY(cb, c)
      DGVS_XOR1(cb, bw, iq)
      /* test `cwb' is the vertices adjacent to `w', limited to `cb' */
      DGVS_AND2(cwb, g->c[w], cb)
      /* `cwb' must be equal to `cb' if `c' is a clique */
      if ( !dgvs_eq(cwb, cb) ) /* not a clique */
        break; /* break the loop prematurely, r != 0 */
      DGVS_XOR1(r, bw, iq)
    }
    if ( !dgvs_nonzero(r) ) { /* if the loop is completed, `c' is a clique */
      if ( dgvs_eq(unvisited, c) ) { /* clique `c' == the rest vertices */
        return ncl;          /* so it is not a separator */
      } else { /* found a clique `c' */
        if (cl != NULL) DGVS_CPY(cl[ncl], c);
        ncl++;
        if (stop1) return 1;
      }
    }
  }
  return ncl;
}


#define dg_cliquesep(g) dg_cliquesep0(g, DGCSEP_DEFAULTMETHOD)

/* test if a graph has a clique separator */
INLINE dgvsref_t dg_cliquesep0(const dg_t *g, int method)
{
  static dgvs_t fs_c[DG_NMAX];
  static dg_t fs[1] = {{0, NULL}}; /* stock fill-in graph */
  static int a[DG_NMAX]; /* a[k] is the kth vertex */
  static dgvs_t cl;
#pragma omp threadprivate(fs_c, fs, a, cl)
  DG_DEFN_(g)

  fs->n = DG_N_;
  fs->c = fs_c;

  method %= DGCSEP_NBASICMETHODS; /* discard higher bits */

  if (method == DGCSEP_MCSP) {
    static int p[DG_NMAX]; /* p[v] is the index of vertex v */
#pragma omp threadprivate(p)
    /* 1. find MCS-P ordering */
    dg_mcsearchp(g, a, p);
    /* 2. compute the fill-in */
    dg_fillin(g, fs, a, p);
    /* 3. clique decomposition (stop after the first clique) */
    return dg_decompcliqueseplow(g, fs, a, &cl, DGCSEP_FIND1) ? cl : 0;
  } else if (method == DGCSEP_MCSM) {
    return dg_minimalordermcsm(g, fs, a, DGCSEP_FIND1, &cl) ? cl : 0;
  } else { /* DGCSEP_LEXM */
    return dg_minimalorderlexm(g, fs, a, DGCSEP_FIND1, &cl) ? cl : 0;
  }
}



/* number of nodes in the clique-separator decomposition */
#define dg_ncsep(g) dg_decompcsep(g, NULL, DGCSEP_DEFAULTMETHOD)

INLINE int dg_decompcsep(const dg_t *g, dgvs_t * RESTRICT cl, int method)
{
  static dgvs_t fs_c[DG_NMAX];
  static dg_t fs[1] = {{0, NULL}}; /* stock fill-in graph */
  static int a[DG_NMAX]; /* a[k] is the kth vertex */
#pragma omp threadprivate(fs_c, fs, a)
  DG_DEFN_(g)

  fs->n = DG_N_;
  fs->c = fs_c;

  method %= DGCSEP_NBASICMETHODS; /* discard higher bits */

  if (method == DGCSEP_MCSM)
    return dg_minimalordermcsm(g, fs, a, DGCSEP_FINDALL, cl);
  else
    return dg_minimalorderlexm(g, fs, a, DGCSEP_FINDALL, cl);
}



#if DGVS_ONEWORD



/* find clique separator of two vertices */
INLINE dgvsref_t dg_csep2(const dg_t *g)
{
  int i;
  dgvs_t ci, maski;
  dgword_t bi, bj;
  DG_DEFN_(g)
  DG_DEFMASKN_()

  /* loop for the first vertex */
  for (i = 1; i < DG_N_; i++) {
    bi = MKBIT(i);
    maski = DG_MASKN_ ^ bi;

    /* loop over vertices connected to `i', and with indices less than `i' */
    for (ci = g->c[i] & MKBITSMASK(i); ci; ci ^= bj) {
      bj = ci & (-ci);
      /* if the cluster is not connected, we have found a clique separator */
      if ( !dg_connectedvs(g, maski ^ bj) )
        return bi ^ bj; /* bi | bj */
    }
  }
  return 0;
}



#else /* !DGVS_ONEWORD */



/* find clique separator of two vertices */
INLINE dgvsref_t dg_csep2(const dg_t *g)
{
  int i;
  dgvs_t ci, maski, maskci;
  dgword_t bi, bj;
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(jq)
  DG_DEFN_(g)
  DG_DEFMASKN_()

  DGVS_CPY(maski, DG_MASKN_)
  DGVS_CLEAR(maskci)
  maskci[0] = 1;

  /* loop for the first vertex */
  for (i = 1; i < DG_N_; i++) {
    DGVS_MKBIT(i, bi, iq) /* bi = MKBIT(i); */
    /* the next two statements; maski = DG_MASKN_ ^ bi; */
    DGVS_XOR1(maski, bi, iq)
    /* maskci: vertices connected to `i', and with indices less than `i' */
    DGVS_AND2(ci, g->c[i], maskci)
    DGVS_OR1(maskci, bi, iq) /* update maskci */

    while ( dgvs_nonzero(ci) ) {
      DGVS_FIRSTBIT(ci, bj, jq)
      DGVS_XOR1(ci, bj, jq) /* remove `bj' from `ci' */
      DGVS_XOR1(maski, bj, jq) /* `maski' now lacks the `j' bit */

      /* if the cluster is not connected, we have found a clique separator */
      if ( !dg_connectedvs(g, maski) ) {
        static dgvs_t csep;
        #pragma omp threadprivate(csep)
        DGVS_CLEAR(csep)
        DGVS_OR1(csep, bi, iq)
        DGVS_OR1(csep, bj, jq)
        return csep;
      }
      DGVS_XOR1(maski, bj, jq) /* add back the `j' bit */
    }
    DGVS_XOR1(maski, bi, iq) /* add back the `i' bit */
  }
  return 0;
}


#endif /* DGVS_ONEWORD */



#if 0
#define dg_csep(g) dg_cseplow(g, 0)

/* test if a graph has a clique separator
 * optimized version of dg_cliquesep(g) */
INLINE dgvsref_t dg_cseplow(const dg_t *g, int ned)
{
  dgvsref_t cc;
  DG_DEFN_(g)

  /* it appears that it is advantageous to test 2-clique separators
   * for loosely-connected diagrams */
  if (ned <= 0) ned = dg_nedges(g);
  if (ned <= DG_N_ * 2 + 3) {
    cc = dg_csep2(g);
    if (cc != 0) return cc;
  }
  return dg_cliquesep(g);
}
#else


#define dg_csep(g) dg_csep0(g, DGCSEP_DEFAULTMETHOD)

/* test if a graph has a clique separator
 * optimized version of dg_cliquesep(g) */
INLINE dgword_t dg_csep0(const dg_t *g, int method)
{
  if (method / DGCSEP_NBASICMETHODS == 0) {
    /* it seems that testing 2-vertex clique separator is
     * open profitable in real sampling, although it is
     * not so for a random graph, as those in testcsep() */
    dgword_t cc = dg_csep2(g);
    if (cc != 0) return cc;
  }
  return dg_cliquesep0(g, method);
}


#endif



#ifdef DGMAP_EXISTS

/* this array is shared due to its large size */
static char *dgcsep_ncl_[DGMAP_NMAX + 1];


/* compute the number of nodes the clique-separator decomposition */
INLINE int dg_ncsep_lookuplow(const dg_t *g, dgword_t c)
{
  int ncs;
  DG_DEFN_(g)

  /* initialize the lookup table */
  if (dgcsep_ncl_[DG_N_] == NULL) {
#pragma omp critical
    {
      /* we test the pointer again because another thread
       * might have allocated the space now */
      if (dgcsep_ncl_[DG_N_] == NULL) {
        int ipr, npr = 1u << (DG_N_ * (DG_N_ - 1) / 2);
        xnew(dgcsep_ncl_[DG_N_], npr);
        for (ipr = 0; ipr < npr; ipr++)
          dgcsep_ncl_[DG_N_][ipr] = (char) (-1);
      }
    }
  }

  ncs = dgcsep_ncl_[DG_N_][c];
  if (ncs < 0) {
    ncs = dg_ncsep(g);
#pragma omp critical
    {
      dgcsep_ncl_[DG_N_][c] = (char) ncs; /* save the value */
    }
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



INLINE void dgcsep_free(void)
{
#ifdef DGMAP_EXISTS
#pragma omp critical
  {
    int i;
    for (i = 0; i <= DGMAP_NMAX; i++)
      if (dgcsep_ncl_[i] != NULL) {
        free(dgcsep_ncl_[i]);
        dgcsep_ncl_[i] = NULL;
      }
  }
#endif /* defined(DGMAP_EXISTS) */
}



#endif /* DGCSEP_H__ */

