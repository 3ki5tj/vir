#ifndef DGRING_H__
#define DGRING_H__
#include "dgrjw.h"



/* compute the number of ring subgraphs
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dgring_spec0(const dg_t *g,
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
       * clique separator
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
    *fb = DG_SC2FB(dgsc_rhiter(DG_N_, 2, 1), ned0);
    return 0;
  }
  return -1;
}



#define dgring_do(g) dgring_iter(g)


/* return the number of ring subgraphs */
INLINE double dgring_iter(const dg_t *g)
{
  DG_DEFN_(g)
  int st[DG_NMAX], top, sttop, root = 0;
  dgvs_t unused, ccp, ccp0, croot, c, masktop, vstmp;
  dgword_t bi;
  DGVS_DEFIQ_(iq)
  double cnt = 0;

  if (DG_N_ <= 2) {
    if (DG_N_ <= 1) return 1;
    else return (double) (int) (DGVS_FIRSTWORD(g->c[1]) & 0x1);
  }
  st[0] = root;
  st[1] = 0;
  top = 1;
  /* unused = DG_MASKN_ ^ MKBIT(root); */
  DGVS_MKBITSMASK(unused, DG_N_)
  DGVS_REMOVE(unused, root)
  DGVS_CPY(croot, g->c[root]) /* vertices adjacent to the `root' */

  DGVS_AND2(ccp, unused, croot); /* ccp: set of unused vertices adjacent to st[top-1] */
  sttop = st[top];

  while (1) {
    DGVS_CPY(ccp0, ccp) /* backup ccp in case pushing failed */
    /* construct a set of vertices that satisfy
     * 1. in ccp: unused && connected to st[top-1]
     * 2. indices > sttop */
    /* c = ccp & MKINVBITSMASK(sttop + 1); */
    DGVS_MKBITSMASK(masktop, sttop + 1)
    DGVS_MINUS2(c, ccp, masktop) /* exclude vertices in masktop */
    if ( dgvs_nonzero(c) ) {
      DGVS_FIRSTLOW(sttop, c, bi, iq) /* choose `sttop' (or bit `bi') from `c' */
      DGVS_XOR1(unused, bi, iq) /* remove b from the unused vertex set */
      DGVS_AND2(ccp, unused, g->c[sttop]) /* ccp = unused & g->c[sttop]; */
      if (ccp) { /* there are still unused vertices */
        if (top < DG_N_ - 2) { /* push */
          st[top++] = sttop; /* save the current top */
          sttop = 0; /* clear the next level */
          continue;
        }
        /* stay on the highest level, fall through
         * ccp should be a single bit representing the last the vertex */
        /* check ccp and the root are adjacent */
        cnt += !dgvs_nonzero( dgvs_and2(vstmp, croot, ccp) );
        // cnt += ((croot & ccp) != 0);
      }
      /* stay on this level */
      DGVS_CPY(ccp, ccp0) /* recover the old ccp */
      DGVS_XOR1(unused, bi, iq) /* add b back to the unused vertices */
      //unused ^= b;
    } else { /* exhausted this level, pop  */
      if (--top == 0) break;
      sttop = st[top];
      DGVS_FLIP(unused, sttop) /* add `b' back to the unused vertices */
      // unused ^= MKBIT(sttop);
      DGVS_AND2(ccp, unused, g->c[st[top - 1]]); /* reconstruct ccp */
    }
  }
  return cnt / 2; /* clockwise & counterclockwise */
}



#define dg_ring(g) dg_ring0(g, NULL, NULL)

/* return the number of ring subgraphs
 * use various techniques to accelerate the calculation
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_ring0(const dg_t *g, int *ned, int *degs)
{
  int err;
  double nr;

  nr = dgring_spec0(g, ned, degs, &err);
  if (err == 0) return nr;
  else return dgring_do(g);
}



/* free all stock objects */
INLINE void dgring_free(void) { }



#endif /* DGRING_H__ */

