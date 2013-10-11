#ifndef DGRING_H__
#define DGRING_H__
#include "dg.h"



/* compute the number of ring subgraphs
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_nring_spec0(const dg_t *g,
    int *ned, int *degs, int *err)
{
  int i, j, n = g->n, ned0, ned1;
  TSTATIC int ldegs[DG_NMAX]; /* local buffer of the degree sequence */
  double x;

  *err = 0;
  /* compute the degrees of all vertices */
  if (degs == NULL) degs = ldegs;
  if (ned == NULL || *ned <= 0) {
    ned0 = dg_degs(g, degs);
    if (ned) *ned = ned0;
  } else {
    ned0 = *ned;
  }

  if (ned0 <= n + 1) {
    if (ned0 == n) return 1;

    /* ned0 == n + 1
     * (a) if the two deg-3 vertices are mutually connected,
     * then there is a ring
     * (b) if the two deg-3 vertices are connected to deg-2 vertices
     * then no subgraph is ring */
    for (i = 0; i < n; i++)
      if (degs[i] == 3) break;
    for (j = i + 1; j < n; j++)
      if (degs[j] == 3) break;
    /* if j >= n, diagram is not biconnected */
    return j < n && dg_linked(g, i, j);
  }

  ned1 = n * (n - 1) / 2 - ned0;
  if (ned1 == 0) { /* full diagram return n!/(2*n) */
    for (x = 1, i = 3; i < n; i++) x *= i;
    return x;
  }

  *err = 1; /* failed */
  return 0;
}



/* return the number of ring subgraphs */
INLINE double dg_nring_direct(const dg_t *g)
{
  int st[DG_NMAX], top, i, n = g->n, root = 0;
  code_t unused, c, b, mask = mkbitsmask(n);
  double cnt = 0;

  if (n <= 2) {
    if (n <= 1) return 1;
    else return (g->c[1] & 0x1) ? 1 : 0;
  }
  st[0] = root;
  st[1] = 0;
  top = 1;
  unused = mask ^ MKBIT(root);
  while (1) {
    /* construct a set of vertices that satisfy
     * 1. unused
     * 2. connected to st[top-1]
     * 3. indices > st[top] */
    if ( (c = unused & g->c[ st[top - 1] ]
            & ~mkbitsmask(st[top] + 1)) != 0 ) {
      i = bitfirstlow(c, &b);
      unused ^= b;
      st[top] = i;
      if (unused) { /* there are still unused vertices */
        st[++top] = 0; /* push, clear the next level */
      } else { /* no vertex left, stay in the highest level */
        if (g->c[root] & b) cnt += 1;
        unused ^= b;
      }
    } else { /* exhausted this level */
      if (--top <= 0) break;
      unused ^= MKBIT(st[top]);
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



#define dg_nring_lookup(g) dg_nring_lookuplow(g->n, dg_getmapid(g))

/* compute the number of ring subgraphs by a look up table */
INLINE double dg_nring_lookuplow(int n, unqid_t id)
{
  static double *nr[DGMAP_NMAX + 1]; /* nr of unique diagrams */

  if (n <= 1) return 1;

#pragma omp critical
  {
    if (nr[n] == NULL) { /* initialize the look-up table */
      dg_t *g;
      dgmap_t *m = dgmap_ + n;
      int k, cnt = 0, nz = 0;
      double *nrn;
      clock_t t0 = clock();

      dgmap_init(m, n);
      xnew(nrn, m->ng);

      /* loop over unique diagrams */
      g = dg_open(n);
      for (cnt = 0, k = 0; k < m->ng; k++) {
        dg_decode(g, &m->first[k]);
        if ( dg_biconnected(g) ) {
          nrn[k] = dg_nring_mixed(g);
          cnt++;
          nz++;
        } else nrn[k] = 0;
      }
      dg_close(g);
      printf("%4d: n %d, computed # of subrings of %d/%d biconnected diagrams, %gs\n",
          inode, n, cnt, nz, 1.*(clock() - t0)/CLOCKS_PER_SEC);
      nr[n] = nrn;
    }
  } /* omp critical */
  return nr[ n ][ id ];
}



#define dg_nring(g) dg_nring0(g, NULL, NULL)

/* return the number of ring subgraphs
 * using various techniques to accelerate the calculation */
INLINE double dg_nring0(const dg_t *g, int *ned, int *degs)
{
  if (g->n <= DGMAP_NMAX) {
    return dg_nring_lookup(g);
  } else {
    return dg_nring_mixed0(g, ned, degs);
  }
}


#endif

