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
INLINE double dg_rhiter(int n, int n0, double sc)
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



/* minimize the diagram of `g' by removing fully-connected vertices */
INLINE dg_t *dg_mintop(dg_t *g)
{
  int i;
  dgword_t mask = MKBITSMASK(g->n);

  /* remove fully-connected vertices until it is no longer biconnected */
  for (i = 0; i < g->n; i++) {
    /* construct a vertex set without i */
    dgword_t vs = mask ^ MKBIT(i);
    /* see if removing a fully-connected vertex leaves
     * the diagram biconnected */
    if ( g->c[i] == vs && dg_biconnectedvs(g, vs) ) {
      dg_remove1(g, g, i);
    } /* otherwise keep searching */
  }
  return g;
}



/* check if a graph is biconnected under the assumption that
 * only vertices in avs can be articulation points */
INLINE int dg_biconnectedavs(const dg_t *g, dgword_t avs)
{
  DG_DEFN_(g);
  DG_DEFMASKN_();
  dgword_t b;

  for (b = 1; b & DG_MASKN_; b <<= 1)
    if ( (b & avs) && !dg_connectedvs(g, DG_MASKN_ ^ b) )
        return 0;
  return 1;
}



#if 0 /* same as the above but slightly slower version */
INLINE int dg_biconnectedavsb(const dg_t *g, dgword_t avs)
{
  DG_DEFN_(g);
  DG_DEFMASKN_();
  dgword_t b;

  for (; avs; avs ^= b) /* `b' is the first nonzero bit */
    if ( !dg_connectedvs(g, DG_MASKN_ ^ (b = avs & (-avs))) )
      return 0;
  return 1;
}
#endif



/* recursively find the star content
 * starting from the edge (i, j + 1) */
INLINE double dg_rhsc_recur(dg_t *g, int sgn, int i, int j)
{
  DG_DEFN_(g);
  DG_DEFMASKN_();
  dgword_t avs, avsi;
  double sc = 0;

  /* find the pair after (i, j) with i < j */
  if (++j >= DG_N_) j = (++i) + 1;
  /* loop over the first vertex i */
  for (; i < DG_N_ - 1; j = (++i) + 1) {
    /* if the degree <= 2, removing an edge connecting i
     * makes the biconnectivity impossible */
    if (dg_deg(g, i) <= 2) continue;
    avsi = DG_MASKN_ ^ MKBIT(i);
    /* loop over the second vertex j */
    for (; j < DG_N_; j++) {
      /* try to remove the edge i-j */
      if (!dg_linked(g, i, j) || dg_deg(g, j) <= 2)
        continue;
      DG_UNLINK(g, i, j);
      /* It is certain that neither i or j is an aritculation points (*)
       * we will only test the connectivity of g without other vertices
       * Prove (*):
       * Since `g' with (i, j) is biconnected, g\{i} is connected
       * so g'\{i} where g' is g \ {(i, j)} is also connected,
       * hence i is not an articulation point of g'  */
      avs = avsi ^ MKBIT(j);
      if ( dg_biconnectedavs(g, avs) ) {
        sc += sgn /* add the diagram without (i, j) */
            + dg_rhsc_recur(g, -sgn, i, j); /* diagrams without (i, j) */
      }
      DG_LINK(g, i, j); /* link back, find subdiagrams with (i, j) */
    }
  }
  return sc;
}



#define dg_rhsc_spec(g, err) dg_rhsc_spec0(g, 0, 1, NULL, NULL, err)

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
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_rhsc_spec0(const dg_t *g, int nocsep, int testcsep,
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
    /* if j >= n, diagram is not biconnected */
    return j < n && !dg_linked(g, i, j);
  }

  /* densely-connected diagrams */
  ned1 = n * (n - 1) / 2 - ned0;
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
  if (testcsep) {
    /* if nocsep, we know for sure there is no clique separator
     * then the hard calculation must be done */
    *err = (nocsep || !dg_cliquesep(g));
  } else { /* if we don't want to test clique separator */
    *err = 1; /* simply show failure */
  }
  return 0;
}



/* compute the star content (SC) of a Ree-Hoover diagram
 * unconnected edge is treated as a wiggly line
 * SC = # of biconnected subgraphs with even edges removed
 *    - # of biconnected subgraphs with odd edges removed */
INLINE double dg_rhsc_directlow(const dg_t *g)
{
  double sc;
  dg_t *g0 = NULL;

  /* first find the minimal set of vertices that
   * contain the wiggly lines */
  g0 = dg_clone(g);
#ifndef N /* mintop is disabled for fixed-size graphs */
  dg_mintop(g0);
#endif
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
INLINE double dg_rhsc_direct0(const dg_t *g, int nocsep, int *ned, int *degs)
{
  double sc;
  int err;

  /* detect special cases when possible */
  sc = dg_rhsc_spec0(g, nocsep, 1, ned, degs, &err);
  if (err == 0) return sc;
  else return dg_rhsc_directlow(g);
}



#ifdef DGMAP_EXISTS
/* compute the star content by a look up table */
INLINE double dg_rhsc_lookup(const dg_t *g)
{
  static double *sc[DGMAP_NMAX + 1];
#pragma omp threadprivate(sc)
  DG_DEFN_(g);
  dgmap_t *m = dgmap_ + DG_N_;
  dgword_t c;

  if (DG_N_ <= 1) return 1;

  if (sc[DG_N_] == NULL) {
    dg_t *g1;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock(), t1;

    if (DG_N_ >= 8) fprintf(stderr, "n %d: initializing...\n", DG_N_);
    dgmap_init(m, DG_N_); /* compute the permutation mapping */
    xnew(sc[DG_N_], m->ng);

    t1 = clock();
    if (DG_N_ >= 8) fprintf(stderr, "n %d: diagram-map initialized %gs\n",
        DG_N_, 1.*(t1 - t0)/CLOCKS_PER_SEC);
    /* loop over unique diagrams */
    g1 = dg_open(DG_N_);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g1, &m->first[k]);
      if ( dg_biconnected(g1) ) {
        sc[DG_N_][k] = dg_rhsc_direct(g1);
        cnt++;
        nz += (fabs(sc[DG_N_][k]) > 0.5);
      } else sc[DG_N_][k] = 0;
    }
    dg_close(g1);
    fprintf(stderr, "n %d, computed star contents of %d/%d biconnected diagrams, %gs\n",
        DG_N_, cnt, nz, 1.*(clock() - t1)/CLOCKS_PER_SEC);
  }
  dg_encode(g, &c);
  return sc[DG_N_][ m->map[c] ]; /* m->map[c] is the id of the unique diagram */
}
#endif /* defined(DGMAP_EXISTS) */


#define dg_rhsc(g) dg_rhsc0(g, 0, NULL, NULL)

/* compute the star content (SC) of a Ree-Hoover diagram
 * see the comments of dg_rhsc_direct0() for details */
INLINE double dg_rhsc0(const dg_t *g, int nocsep, int *ned, int *degs)
{
  double sc;
  int err;

#ifdef DGMAP_EXISTS
  if (g->n <= DGMAP_NMAX)
    return dg_rhsc_lookup(g);
#endif /* defined(DGMAP_EXISTS) */

  sc = dg_rhsc_spec0(g, nocsep, 1, ned, degs, &err);
  if (err == 0) return sc;
  else return dg_rhsc_directlow(g);
}



#endif

