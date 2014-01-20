#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content by the direct method */



#include <time.h>
#include <limits.h>
#include "dgcsep.h"
#include "dgutil.h"



/* conversion between star-content and hard-sphere weight */
#define DG_SC2FB(sc, ned) ((sc) * (1 - (ned % 2) * 2))
#define DG_FB2SC(fb, ned) DG_SC2FB(fb, ned)



/* Ree-Hoover formula for the star content of a larger diagram of
 * n vertices from that of a smaller diagram of n-1 vertices */
INLINE double dgsc_rhiter(int n, int n0, double sc)
{
  int i;

  if (n < n0) return 0;
  for (i = n0; i < n; i++) /* Ree and Hoover formula */
    sc *= (i - 1) * (i % 2 ? -1 : 1);
  /* e.g. if SC(2) = 1, then
   * n:   2  3  4  5  6   7    8     9    10
   * SC:  1  1 -2 -6 24 120 -720 -5040 40320 */
  return sc;
}



#define dgsc_spec(g, err) \
  dgsc_spec0(g, DGCSEP_DEFAULTMETHOD, NULL, NULL, err)

/* compute the star content (SC) of a Ree-Hoover diagram in special cases
 * only use cheap strategies to deduce the SC
 * assuming the diagram is biconnected
 * if successful, *err = 0, otherwise *err = 1
 *   if *err = 0, the diagram has a clique separator if the returned
 *                value of SC is zero
 *   if *err = 1, the diagram has no clique separator on return
 * csepmethod = 0 means not to test clique separators
 *            = 1 means to test clique separators thoroughly
 *            = 2 means to test clique separators by maximal cardinality search
 *              which can produce false negatives
 * *ned: number of edges; degs: unsorted degree sequence
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dgsc_spec0(const dg_t *g, int csepmethod,
    int *ned, int *degs, int *err)
{
  int i, j, n = g->n, ned0, ned1;
  static int ldegs[DG_NMAX]; /* local buffer for the degree sequence */
#pragma omp threadprivate(ldegs)

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
    if (ned0 == n) /* ring diagram or intertwined rings (see below) */
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
    /* if j >= n, diagram is not biconnected */
    return j < n && !dg_linked(g, i, j);
  }

  /* densely-connected diagrams */
  ned1 = n * (n - 1) / 2 - ned0;
  if (ned1 < 3) {
    if (ned1 == 0) { /* fully connected */
      return dgsc_rhiter(n, 2, 1);
    } else if (ned1 == 1) { /* always has a clique separator */
      return 0;
    } else { /* ned1 == 2 */
      /* if there are two wiggly lines, the SC is nonzero only if
       * the two wiggly lines are not connected to the same vertex */
      for (i = 0; i < n; i++)
        if (degs[i] <= n - 3) return 0;
      return dgsc_rhiter(n, 4, 1);
    }
  }

  /* general case: try to find a clique separator */
  if (csepmethod) {
    /* csepmethod can be 1 or 2, in the latter case
     * we use dg_csep(g, 2), which invokes the maximal cardinality
     * search; this may yield false negative result */
    /* if there is a clique separator, dg_csep0() returns nonzero
     * and we know fb == 0, and there is no error */
    *err = !dg_csep0(g, csepmethod);
  } else { /* if we don't want to test clique separator */
    *err = 1; /* simply show failure */
  }
  return 0;
}




/* minimize the diagram of `g' by removing fully-connected vertices */
INLINE dg_t *dgsc_mintop(dg_t *g)
{
  int i;
  dgvs_t vs;
  dgword_t bi;
  DGVS_DEFIQ_(iq)

  /* remove fully-connected vertices until it is no longer biconnected
   * note, g->n changes in the loop */
  for (i = 0; i < g->n; i++) {
    /* construct a vertex set without i */
    DGVS_MKBIT(i, bi, iq)
    DGVS_MKBITSMASK(vs, g->n);
    DGVS_XOR1(vs, bi, iq) /* remove `bi' from `vs' */
    /* see if removing a fully-connected vertex leaves
     * the diagram biconnected */
    if ( dgvs_eq(g->c[i], vs) && dg_biconnectedvs(g, vs) ) {
      //dg_print(g);
      //printf("removing %d/%d\n", i, g->n);
      //dgvs_printn(vs, "vs");
      dg_remove1(g, g, i);
      //dg_print(g);
    } /* otherwise keep searching */
    DGVS_XOR1(vs, bi, iq) /* add `bi' back to `vs' */
  }
  return g;
}



/* recursively find the star content
 * starting from the edge (i, j + 1) */
INLINE double dgsc_recur(dg_t *g, int sgn, int i, int j)
{
  //dgvs_t avs;
  dgword_t bi, bj;
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(jq)
  double sc = 0;
  DG_DEFN_(g)

  /* find the pair after (i, j) with i < j */
  if (++j >= DG_N_) j = (++i) + 1;
  //DGVS_MKBITSMASK(avs, DG_N_)
  //dg_print(g); printf("n %d\n", DG_N_); getchar();
  /* loop over the first vertex i */
  for (; i < DG_N_ - 1; j = (++i) + 1) {
    /* if the degree <= 2, removing an edge connecting i
     * makes the biconnectivity impossible */
    if (dg_deg(g, i) <= 2) continue;
    DGVS_MKBIT(i, bi, iq)
    //DGVS_XOR1(avs, bi, iq) /* remove `bi' from `avs' */
    /* loop over the second vertex j */
    for (; j < DG_N_; j++) {
      /* try to remove the edge i-j */
      if ( !DGVS_HASBIT(g->c[j], bi, iq) )
        continue;
      if (dg_deg(g, j) <= 2)
        continue;
      //printf("n %d, before removing (%d, %d), bi %x bj %x ci %x, cj %x\n", DG_N_, i, j, bi, bj, g->c[i], g->c[j]);
      DGVS_MINUS1(g->c[j], bi, iq) /* remove i from c[j] */
      //g->c[j] &= ~bi;
      DGVS_MKBIT(j, bj, jq)
      DGVS_MINUS1(g->c[i], bj, jq) /* remove j from c[i] */
      //dg_print(g);
      //printf("n %d, removing (%d, %d), bi %x bj %x ci %x, cj %x\n", DG_N_, i, j, bi, bj, g->c[i], g->c[j]); getchar();
      //g->c[i] &= ~bj;
      /* It is certain that neither i or j is an articulation points (*)
       * we will only test the connectivity of g without other vertices
       * Prove (*):
       * Since `g' with (i, j) is biconnected, g\{i} is connected
       * so g'\{i} where g' is g \ {(i, j)} is also connected,
       * hence i is not an articulation point of g'  */
      //DGVS_XOR1(avs, bj, jq) /* remove `bj' from `avs' */
      //if ( dg_biconnectedavs(g, avs) ) {
      if ( dg_biconnected(g) ) {
        sc += sgn /* add the diagram without (i, j) */
            + dgsc_recur(g, -sgn, i, j); /* diagrams without (i, j) */
      }
      //DGVS_XOR1(avs, bj, jq) /* add back `bj' to `avs' */
      /* link back, find subdiagrams with (i, j) */
      DGVS_OR1(g->c[i], bj, jq)
      DGVS_OR1(g->c[j], bi, iq)
      //g->c[i] |= bj;
      //g->c[j] |= bi;
    }
    //DGVS_XOR1(avs, bi, iq) /* add `bi' back to `avs' */
  }
  return sc;
}



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed
 *    - # of biconnected subgraphs with odd edges removed */
INLINE double dgsc_do(const dg_t *g)
{
  double sc;
  static dgvs_t g0_c[DG_NMAX];
  static dg_t g0[1] = {{0, NULL}};
#pragma omp threadprivate(g0_c, g0)

  /* first find the minimal set of vertices that
   * contain the wiggly lines */
  g0->n = g->n;
  g0->c = g0_c;
  dg_copy(g0, g);
#ifndef N /* mintop is disabled for fixed-size graphs */
  dgsc_mintop(g0);
#endif
  sc = 1 + dgsc_recur(g0, -1, 0, 0);
  /* use the Ree-Hoover formula to go from n0 to n */
  sc = dgsc_rhiter(g->n, g0->n, sc);
  return sc;
}



#define dg_sc(g) dg_sc0(g, DGCSEP_DEFAULTMETHOD, NULL, NULL)

/* directly compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed -
 *      # of biconnected subgraphs with odd edges removed
 * nocsep = 1 means if the graph has been tested with no clique separator
 *        = 0 means it MAY have clique separators
 * *ned: number of edges; degs: degree sequence */
INLINE double dg_sc0(const dg_t *g, int csepmethod,
    int *ned, int *degs)
{
  double sc;
  int err;

  /* detect special cases when possible */
  sc = dgsc_spec0(g, csepmethod, ned, degs, &err);
  if (err == 0) return sc;
  else return dgsc_do(g);
}



#define dg_rhsc(g) dg_rhsc0(g, DGCSEP_DEFAULTMETHOD, NULL, NULL)

/* compute the star content (SC) of a Ree-Hoover diagram
 * see the comments of dg_rhsc_direct0() for details */
INLINE double dg_rhsc0(const dg_t *g, int csepmethod,
    int *ned, int *degs)
{
  double sc;
  int err;

  sc = dgsc_spec0(g, csepmethod, ned, degs, &err);
  if (err == 0) return sc;
  else return dgsc_do(g);
}



#endif /* !defined(DGSC_H__) */

