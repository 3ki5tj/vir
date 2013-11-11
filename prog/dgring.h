#ifndef DGRING_H__
#define DGRING_H__
#include "dgrjw.h"



/* compute the number of ring subgraphs
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_nring_spec0(const dg_t *g,
    int *ned, int *degs, int *err)
{
  int i, j, ned0, ned1;
  DG_DEFN_(g)
  double x;
  static int ldegs[DG_NMAX]; /* local buffer of the degree sequence */
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

  if (ned0 <= DG_N_ + 1) {
    if (ned0 == DG_N_) return 1;

    /* ned0 == n + 1
     * (a) if the two deg-3 vertices are mutually connected,
     * then there is a ring
     * (b) if the two deg-3 vertices are connected to deg-2 vertices
     * then no subgraph is ring */
    for (i = 0; i < DG_N_; i++)
      if (degs[i] == 3) break;
    for (j = i + 1; j < DG_N_; j++)
      if (degs[j] == 3) break;
    /* if j >= n, diagram is not biconnected */
    return j < DG_N_ && dg_linked(g, i, j);
  }

  ned1 = DG_N_ * (DG_N_ - 1) / 2 - ned0;
  if (ned1 == 0) { /* full diagram return n!/(2*n) */
    for (x = 1, i = 3; i < DG_N_; i++) x *= i;
    return x;
  }

  *err = 1; /* failed */
  return 0;
}



#define dg_fbnr_spec(g, fb, nr) dg_fbnr_spec0(g, fb, nr, NULL, NULL)

/* compute fb and nr in special cases
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE int dg_fbnr_spec0(const dg_t *g, double *fb, double *nr,
    int *ned, int *degs)
{
  int i, j, ned0, ned1;
  DG_DEFN_(g)
  double x;
  static int ldegs[DG_NMAX]; /* local buffer of the degree sequence */
#pragma omp threadprivate(ldegs)

  *fb = *nr = 0;
  /* compute the degrees of all vertices */
  if (degs == NULL) degs = ldegs;
  if (ned == NULL || *ned <= 0) {
    ned0 = dg_degs(g, degs);
    if (ned) *ned = ned0;
  } else {
    ned0 = *ned;
  }

  if (ned0 <= DG_N_ + 1) {
    if (ned0 == DG_N_) {
      *fb = DG_SC2FB(1, ned0);
      *nr = 1;
    } else {
      /* ned0 == n + 1
       * (a) if the two deg-3 vertices are mutually connected,
       * then there is a ring, but the two deg-3 vertices is a
       * clique separtor
       * (b) if the two deg-3 vertices are connected to deg-2 vertices
       * then no subgraph is ring, and sc = 1 */
      for (i = 0; i < DG_N_; i++)
        if (degs[i] == 3) break;
      for (j = i + 1; j < DG_N_; j++)
        if (degs[j] == 3) break;
      /* if j >= n, diagram is not biconnected */
      if (j >= DG_N_) return -1;
      if (dg_linked(g, i, j)) { /* case (a) */
        *fb = 0;
        *nr = 1;
      } else { /* case (b) */
        *fb = DG_SC2FB(1, ned0);
        *nr = 0;
      }
    }
    return 0;
  }

  ned1 = DG_N_ * (DG_N_ - 1) / 2 - ned0;
  if (ned1 == 0) { /* full diagram nr = n!/(2*n), fb */
    for (x = 1, i = 3; i < DG_N_; i++) x *= i;
    *nr = x;
    *fb = DG_SC2FB(dg_rhiter(DG_N_, 2, 1), ned0);
    return 0;
  }
  return -1;
}



/* return the number of ring subgraphs */
INLINE double dg_nring_direct(const dg_t *g)
{
  DG_DEFN_(g)
  DG_DEFMASKN_()
  int st[DG_NMAX], top, sttop, root = 0;
  dgword_t unused, c, b, croot, ccp;
  double cnt = 0;

  if (DG_N_ <= 2) {
    if (DG_N_ <= 1) return 1;
    else return (double) (int) (g->c[1] & 0x1);
  }
  st[0] = root;
  st[1] = 0;
  top = 1;
  unused = DG_MASKN_ ^ MKBIT(root);
  croot = g->c[root];

  ccp = unused & croot; /* ccp: set of unused vertices adjacent to st[top-1] */
  sttop = st[top];

  while (1) {
    dgword_t ccp0 = ccp; /* backup ccp in case pushing failed */
    /* construct a set of vertices that satisfy
     * 1. in ccp: unused && connected to st[top-1]
     * 2. indices > sttop */
    c = ccp & MKINVBITSMASK(sttop + 1);
    if ( c != 0 ) {
      BITFIRSTLOW(sttop, c, b);
      unused ^= b; /* remove b from the unused vertex set */
      ccp = unused & g->c[sttop];
      if (ccp) { /* there are still unused vertices */
        if (top < DG_N_ - 2) { /* push */
          st[top++] = sttop; /* save the current top */
          sttop = 0; /* clear the next level */
          continue;
        }
        /* stay on the highest level, fall through
         * ccp should be a single bit representing the last the vertex */
        cnt += ((croot & ccp) != 0); /* check ccp and the root are adjacent */
      }
      /* stay on this level */
      ccp = ccp0; /* recover the old ccp */
      unused ^= b; /* add b back to the unused vertices */
    } else { /* exhausted this level, pop  */
      if (--top == 0) break;
      sttop = st[top];
      unused ^= MKBIT(sttop); /* add b back to the unused vertices */
      ccp = unused & g->c[st[top - 1]]; /* reconstruct ccp */
    }
  }
  return cnt / 2; /* clockwise & counterclockwise */
}



#define dg_nring_mixed(g) dg_nring_mixed0(g, NULL, NULL)

/* return the number of ring subgraphs
 * use various techniques to accelerate the calculation
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_nring_mixed0(const dg_t *g, int *ned, int *degs)
{
  int err;
  double nr;

  nr = dg_nring_spec0(g, ned, degs, &err);
  if (err == 0) return nr;
  else return dg_nring_direct(g);
}



#ifdef DGMAP_EXISTS
#define dg_nring_lookup(g) dg_nring_lookuplow(g->n, dg_getmapid(g))

static double *dgmap_nr_[DGMAP_NMAX + 1]; /* nr of unique diagrams */
#pragma omp threadprivate(dgmap_nr_)

/* compute the number of ring subgraphs by a look up table */
INLINE double dg_nring_lookuplow(int n, unqid_t id)
{
  if (DG_N_ <= 1) return 1;

  /* initialize the look-up table */
  if (dgmap_nr_[DG_N_] == NULL) {
    dg_t *g;
    dgmap_t *m = dgmap_ + DG_N_;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock();

    dgmap_init(m, DG_N_);
    xnew(dgmap_nr_[DG_N_], m->ng);

    /* loop over unique diagrams */
    g = dg_open(DG_N_);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g, &m->first[k]);
      if ( dg_biconnected(g) ) {
        dgmap_nr_[DG_N_][k] = dg_nring_mixed(g);
        cnt++;
        nz++;
      } else dgmap_nr_[DG_N_][k] = 0;
    }
    dg_close(g);
    fprintf(stderr, "%4d: n %d, computed # of subrings of %d/%d biconnected diagrams, %gs\n",
          inode, DG_N_, cnt, nz, 1.*(clock() - t0)/CLOCKS_PER_SEC);
  }
  return dgmap_nr_[ DG_N_ ][ id ];
}
#endif /* defined(DGMAP_EXISTS) */


#define dg_nring(g) dg_nring0(g, NULL, NULL)

/* return the number of ring subgraphs
 * using various techniques to accelerate the calculation */
INLINE double dg_nring0(const dg_t *g, int *ned, int *degs)
{
#ifdef DGMAP_EXISTS
  if (g->n <= DGMAP_NMAX)
    return dg_nring_lookup(g);
#endif /* defined(DGMAP_EXITS) */
  return dg_nring_mixed0(g, ned, degs);
}



/* free all stock objects */
INLINE void dgring_free(void)
{
#ifdef DGMAP_EXISTS
  int k;

  /* dgmap_nr_[] is private */
  for (k = 0; k <= DGMAP_NMAX; k++)
    if (dgmap_nr_[k] != NULL)
      free(dgmap_nr_[k]);
#endif /* defined(DGMAP_EXISTS) */
}



#endif /* DGRING_H__ */

