#ifndef DGSC__
#define DGSC__
/* compute the Ree-Hoover star content by the direct method */



#include <time.h>
#include <limits.h>
#include "dgutil.h"



/* conversion between star-content and hard-sphere weight */
#define DG_SC2FB(sc, ned) ((sc) * (1 - (ned % 2) * 2))
#define DG_FB2SC(fb, ned) DG_SC2FB(fb, ned)



/* Ree-Hoover formula for the signed star content of a larger diagram of
 * n vertices from that of a smaller diagram of n-1 vertices */
INLINE double dgsc_rhiter(int n, int n0, double fb)
{
  int i;

  if (n < n0) return 0;
  for (i = n0; i < n; i++) /* Ree and Hoover formula */
    fb *= i - 1;
  /* the following works  for the star content
     which requires a sign flipping ever two n
      sc *= (i - 1) * (i % 2 ? -1 : 1);
   */
  /* e.g. if SC(2) = 1, then
   * n:     2  3  4  5   6    7    8     9     10
   * SC:    1  1 -2 -6  24  120 -720 -5040  40320
   * ned:     +2 +3 +4  +5   +6   +7    +8     +9
   * SSC:  -1 -1 -2 -6 -24 -120 -720 -5040 -40320
   * */
  return fb;
}



#define dg_fbnr_spec(g, nr, err) dg_fbnr_spec0(g, nr, NULL, NULL, err)

/* compute the signed star content (fb) and ring content (nr)
 * of a Ree-Hoover diagram in special cases
 * only use cheap strategies to deduce the SSC
 * assuming the diagram is biconnected
 * if successful, *err = 0, otherwise *err = 1
 * *ned: number of edges; degs: unsorted degree sequence
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
static double dg_fbnr_spec0(const dg_t *g, double *nr,
    int *ned, int *degs, int *err)
{
  int i, j, n = g->n, ned1;
  double fb, ring;

  *err = 0;
  DG_CALC_DEGS(ned, degs, dg_nedges_, dg_degs_); /* prepare degrees */
  if (nr == NULL) nr = &ring;

  /* loosely connected diagrams */
  if (*ned <= n + 1) {
    if (*ned == n) { /* ring diagram or intertwined rings (see below) */
      *nr = 1;
      fb = DG_SC2FB(1, *ned);
    } else {
      /* *ned == n + 1
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
      if (j < n && !dg_linked(g, i, j)) {
        *nr = 0;
        fb = DG_SC2FB(1, *ned);
      } else {
        *nr = 1;
        fb = 0;
      }
    }
    return fb;
  }

  /* densely-connected diagrams */
  ned1 = n * (n - 1) / 2 - *ned;
  if (ned1 < 3) {
    /* ring content
     * fully connected
     *    n!/(2n) = (n-1)!/2
     * one adjacent wiggly line
     *    n!/(2n) - (n-2)! = (n-2)! (n - 3) / 2
     * where the (n-2)! accounts rings that pass through the missing edge
     * two adjacent wiggly lines
     *    n!/(2n) - (n-2)!*2 + 2 (n - 3)!
     * where the last sum accounts for edges pass through both missing edges
     * they have been subtracted twice in counting edges pass through either edge
     * so we add them back here
     * two nonadjacent wiggly lines:
     *    n!/(2n) - (n-2)!*2 + (n-3)!
     * */
    double *fact = dg_facts_;
    DG_CALC_FACTS()
    if (ned1 == 0) { /* fully connected */
      *nr = fact[n - 1] * .5;
      fb = dgsc_rhiter(n, 2, -1);
    } else if (ned1 == 1) { /* always has a clique separator */
      *nr = fact[n - 2] * (.5*n - 1.5);
      fb = 0;
    } else { /* ned1 == 2 */
      /* if there are two wiggly lines, the SC is nonzero only if
       * the two wiggly lines are not adjacent to the same vertex */
      for (i = 0; i < n; i++)
        if (degs[i] < n - 2) break; /* two wiggly lines are adjacent to i */

      if (i < n) { /* two adjacent wiggly lines */
        *nr = fact[n - 2] * (.5 * n - 2.5) + fact[n - 3];
        fb = 0;
      } else { /* two nonadjacent wiggly lines */
        /* there are two ways of linking the two nonadjacent edges a~b and c~d
         * we may connect a and c by a path of i vertices and
         * connect b and d by a path of n - 4 - i vertices
         * or
         * connect a and d by a path of i vertices and
         * and connect b and c by a path of n - 4 - i vertices */
        /* more generally, for k mutually nonadjacent wiggly lines
         * we have (n - 1)!1F1(-k, -n + 1, -2) / 2 */
        *nr = fact[n - 2] * (.5 * n - 2.5) + 2 * fact[n - 3];
        fb = dgsc_rhiter(n, 4, 1);
      }
    }
    return fb;
  }

  *err = 1; /* failure */
  return 0;
}



/* minimize the diagram of `g' by removing fully-connected vertices */
static dg_t *dgsc_mintop(dg_t *g)
{
  int i;
  dgvs_t vs;
  dgword_t bi;
  DGVS_DEFIQ_(iq)

  /* remove fully-connected vertices until it is no longer biconnected
   * note, g->n changes in the loop */
  for (i = 0; i < g->n; ) {
    /* construct a vertex set without i */
    DGVS_MKBIT(i, bi, iq)
    DGVS_MKBITSMASK(vs, g->n);
    DGVS_XOR1(vs, bi, iq) /* remove `bi' from `vs' */
    /* see if removing a fully-connected vertex leaves
     * the diagram biconnected */
    if ( dgvs_eq(g->c[i], vs) && dg_biconnectedvs(g, vs) ) {
      dg_remove1(g, g, i);
    } else { /* otherwise keep searching */
      i++;
    }
  }
  return g;
}



/* recursively find the star content
 * starting from the edge (i, j + 1) */
static double dgsc_recur(dg_t *g, double *nr, int ned, int i, int j)
{
  static double ring;
#pragma omp threadprivate(ring)
  dgword_t bi, bj;
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(jq)
  double fb = 1 - ned % 2 * 2;
  int n = g->n;

  if (nr == NULL) { nr = &ring; ring = 0; }
  if (i == 0 && j == 0) *nr = 0;
  if (ned == n) *nr += 1;

  /* find the pair after (i, j) with i < j */
  if (++j >= n) j = (++i) + 1;
  /* loop over the first vertex i */
  for (; i < n - 1; j = (++i) + 1) {
    /* if the degree <= 2, removing an edge connecting i
     * makes the biconnectivity impossible */
    if (dg_deg(g, i) < 3) continue;
    DGVS_MKBIT(i, bi, iq)
    /* loop over the second vertex j */
    for (; j < DG_N_; j++) {
      /* try to remove the edge i-j */
      if ( !DGVS_HASBIT(g->c[j], bi, iq) || dg_deg(g, j) < 3 )
        continue;
      DGVS_MINUS1(g->c[j], bi, iq) /* remove i from c[j] */
      DGVS_MKBIT(j, bj, jq)
      DGVS_MINUS1(g->c[i], bj, jq) /* remove j from c[i] */
      /* It is certain that neither i or j is an articulation points (*)
       * we will only test the connectivity of g without other vertices
       * Prove (*):
       * Since `g' with (i, j) is biconnected, g\{i} is connected
       * so g'\{i} where g' is g \ {(i, j)} is also connected,
       * hence i is not an articulation point of g'  */
      if ( dg_biconnected(g) )
        fb += dgsc_recur(g, nr, ned - 1, i, j); /* diagrams without (i, j) */
      /* link back, find subdiagrams with (i, j) */
      DGVS_OR1(g->c[i], bj, jq)
      DGVS_OR1(g->c[j], bi, iq)
    }
  }
  return fb;
}



#define dgsc_iter(g, nr) dgsc_iter0(g, nr, NULL, NULL)

/* iteratively compute the star content,
 * on return g should be the same */
static double dgsc_iter0(dg_t *g, double *nr, int *ned, int *degs)
{
  static int ed[DG_NMAX * (DG_NMAX - 1)/2][2]; /* edges */
  static int st[DG_NMAX * (DG_NMAX - 1)/2 + 2]; /* state */
#pragma omp threadprivate(ed, st)
  int n = g->n, top, ied, med, vi, vj, degi, degj;
  double fb, ring;

  /* collect edges */
  DG_CALC_DEGS(ned, degs, dg_nedges_, dg_degs_)
  for (med = 0, vi = 0; vi < n - 1; vi++)
    if ( degs[vi] > 2 )
      for (vj = vi + 1; vj < n; vj++)
        if ( degs[vj] > 2 && dg_linked(g, vi, vj) )
          ed[med][0] = vi, ed[med][1] = vj, med++;
  if (nr == NULL) nr = &ring;

  st[ top = 0 ] = -1;
  fb = 1 - *ned % 2 * 2;
  *nr = (*ned == n) ? 1 : 0;
  while (top >= 0) {
    /* search over remaining edges */
    for (ied = st[top] + 1; ied < med; ied++) {
      /* removing an edge adjacent to a degree-2 vertex
       * makes a graph not biconnected */
      if ( (degi = degs[ vi = ed[ied][0] ]) < 3
        || (degj = degs[ vj = ed[ied][1] ]) < 3 )
        continue;
      dg_unlink(g, vi, vj);
      if ( dg_biconnected(g) ) { /* push */
        degs[vi] = degi - 1, degs[vj] = degj - 1;
        st[top] = ied;
        --(*ned);
        fb += 1 - *ned % 2 * 2;
        if (*ned == n) *nr += 1;
        break;
      }
      dg_link(g, vi, vj);
    }
    if (ied < med) { /* try to push */
      if (++top < med) {
        /* edges in the stack are arranged in ascending order, so the edge
         * on the next position must have a greater index than this edge */
        st[top] = ied; /* clean up the next level */
        continue;
      } /* otherwise, fall through and pop */
    }
    --top; /* pop */
    if (top >= 0 && (ied = st[top]) >= 0) { /* erase the edge */
      vi = ed[ied][0]; vj = ed[ied][1];
      dg_link(g, vi, vj);
      degs[vi]++, degs[vj]++;
      (*ned)++;
    }
  }
  return fb;
}



enum {
  DGSC_NULLMETHOD = 0,
  DGSC_ITER = 0x100,
  DGSC_RECUR = 0x200,
  DGSC_MASKBASICMETHODS = 0xff00,
  DGSC_NOSPEC = 0x10000 /* don't detect special cases */
};

#ifndef DGSC_DEFAULTMETHOD
#define DGSC_DEFAULTMETHOD DGSC_ITER
#endif

#define dgsc_fb(g) dgsc_fb0(g, DGSC_DEFAULTMETHOD, NULL, NULL)
#define dgsc_fb0(g, method, ned, degs) dgsc_fbnr0(g, NULL, method, ned, degs)
#define dgsc_fbnr(g, nr) dgsc_fbnr0(g, nr, DGSC_DEFAULTMETHOD, NULL, NULL)

/* compute the signed star content (SSC) of a Ree-Hoover diagram
 * if nr != NULL, the ring content is also computed
 * an unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed
 *    - # of biconnected subgraphs with odd edges removed
 * SSC = SC * (-)^(# of edges) */
static double dgsc_fbnr0(const dg_t *g, double *nr, int method,
    int *ned, int *degs)
{
  double fb = 0, ring;
  static dgvs_t g0_c[DG_NMAX];
  static dg_t g0[1] = {{0, NULL}};
#pragma omp threadprivate(g0_c, g0)

  /* first find the minimal set of vertices that
   * contain the wiggly lines */
  g0->n = g->n;
  g0->c = g0_c;
  dg_copy(g0, g);
  /* since we do not have a simple recursion formula for the ring content,
   * the vertex removal technique is disabled when the ring content is needed */
  if (nr == NULL) {
    dgsc_mintop(g0);
    nr = &ring; /* still provide a buffer, but *nr is for g0, not g */
  }
  *nr = 0;

  if ( (method & DGSC_MASKBASICMETHODS) == DGSC_ITER ) {
    fb = dgsc_iter0(g0, nr, ned, degs);
  } else if ( (method & DGSC_MASKBASICMETHODS) == DGSC_RECUR) {
    DG_CALC_DEGS(ned, degs, dg_nedges_, dg_degs_)
    fb = dgsc_recur(g0, nr, *ned, 0, 0);
  }

  /* use the Ree-Hoover formula to go from n0 to n */
  return dgsc_rhiter(g->n, g0->n, fb);
}



#define dg_ssc(g) dg_ssc0(g, 0, NULL, NULL)

/* compute the signed star content (SSC or fb) of a Ree-Hoover diagram
 * same function as dg_fb0() in dgrjw.h, but not as powerful
 * no clique separator is detected and Wheatley's method is unused */
INLINE double dg_ssc0(const dg_t *g, int method, int *ned, int *degs)
{
  double fb;
  int err;

  if ( !(method & DGSC_NOSPEC) ) {
    fb = dg_fbnr_spec0(g, NULL, ned, degs, &err);
    if ( !err ) return fb;
  }
  fb = dgsc_fbnr0(g, NULL, method, ned, degs);
  return fb;
}



#endif /* !defined(DGSC_H__) */

