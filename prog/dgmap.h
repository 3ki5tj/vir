#ifndef DGMAP_H__
#define DGMAP_H__
#include "dg.h"


/* the regular lookup table for n <= DGMAP_NMAX */
#ifndef DGMAP_NMAX
#define DGMAP_NMAX 8
#endif

/* enable the regular lookup table if N is not fixed or it is small enough */
#if !defined(N) || N <= DGMAP_NMAX
  #ifndef DGMAP_EXISTS
  #define DGMAP_EXISTS 1
  #endif
#else
  #ifdef DGMAP_EXISTS
  #undef DGMAP_EXISTS
  #endif
#endif


#ifdef DGMAP_EXISTS
typedef short unqid_t;

typedef struct {
  int ng; /* number of unique diagrams */
  unqid_t *map; /* map a diagram to the unique diagram
                   `short' works for n <= 8 */
  dgword_t *first; /* index of the first unique diagram */
} dgmap_t;



dgmap_t dgmap_[DGMAP_NMAX + 1];



/* compute all permutations of n */
INLINE int dgmap_getperm(int n, int **pp)
{
  int np, npp, ipp, i, top;
  static int st[DGMAP_NMAX + 2], used[DGMAP_NMAX + 2];
#pragma omp threadprivate(st, used)

  for (np = 1, i = 1; i <= n; i++) np *= i;
  npp = np * n; /* each permutation needs n numbers */
  xnew(*pp, npp);
  /* add all permutations of the diagram */
  for (i = 0; i < n; i++) {
    st[i] = -1;
    used[i] = 0;
  }
  for (ipp = 0, top = 0; ; ) {
    if (top >= n) {
      /* found an index permutation, build the graph */
      die_if (ipp >= npp, "ipp %d >= npp %d, n %d\n", ipp, npp, n);
      for (i = 0; i < n; i++) {
        die_if (st[i] < 0 || st[i] >= n, "bad %d/%d %d\n", i, n, st[i]);
        (*pp)[ipp++] = st[i];
      }
    } else {
      /* select the number on level top */
      for (i = st[top]; ++i < n; )
        if ( !used[i] ) break;
      if (i < n) {
        used[i] = 1;
        st[top] = i;
        st[++top] = -1; /* clear the next level */
        continue;
      }
    }
    /* exhausted this level */
    if (--top < 0) break;
    used[ st[top] ] = 0;
  }
  return np;
}


/* compute initial diagram map */
INLINE int dgmap_init(dgmap_t *m, int n)
{
  dgword_t c, c1, ng, *masks, *ms;
  int ipm, npm, *pm, i, j, ipr, npr, gid;
  clock_t t0;
  /* number of unique diagrams */
  static int nfirst[DGMAP_NMAX + 1] = {1,
    1, 2, 4, 11, 34, 156, 1044, 12346};

  die_if (n > DGMAP_NMAX || n <= 0, "bad n %d\n", n);
  /* if any thread sees `m->ng > 0' it is already initialized */
  if (m->ng > 0) return 0; /* already initialized */

#ifdef _OPENMP
#pragma omp critical
  {
    /* the `if' clause is necessary, it avoids multiple threads
     * doing the initialization.  Suppose multiple threads reach here
     * and see m->ng == 0, the first thread enters the block and
     * does the initialization while the others wait.  Then, when the
     * second thread enters the block, it can skip the steps only
     * if it see that `m->ng' has been initialized */
    if (m->ng <= 0) {
#endif /* _OPENMP */
      t0 = clock();
      if (n >= 8) fprintf(stderr, "%4d: n %d: initializing the diagram map\n", inode, n);

      npr = n * (n - 1) / 2;
      ng = (int)( 1u << npr );
      xnew(m->map, ng);
      fprintf(stderr, "%4d: dgmap allocated %gMB for n %d\n",
          inode, (sdgword_t) ng * sizeof(m->map[0]) / (1024.*1024), n);
      for (c = 0; c < ng; c++) m->map[c] = -1;
      xnew(m->first, nfirst[n]);

      if (n == 1) {
        m->map[0] = 0;
        m->first[0] = 0;
        goto END;
      }

      /* compute all permutations */
      npm = dgmap_getperm(DG_N_, &pm);
      /* for each permutation compute the mask of each vertex pair */
      xnew(masks, npm * npr);
      for (ipm = 0; ipm < npm; ipm++)
        for (ipr = 0, i = 0; i < DG_N_ - 1; i++)
          for (j = i + 1; j < DG_N_; j++, ipr++)
            masks[ipm * npr + ipr] /* code bit of the pair (i, j) */
              = MKBIT( getpairindex(
                    pm[ipm*DG_N_ + i], pm[ipm*DG_N_ + j], DG_N_
                    ) );
      free(pm);

      /* loop over all diagrams */
      for (gid = 0, c = 0; c < ng; c++) {
        if (m->map[c] >= 0) continue;
        die_if (gid >= nfirst[n], "too many diagrams %d\n", gid);
        m->first[gid] = c;
        /* add all permutations of the diagram */
        for (ms = masks, ipm = 0; ipm < npm; ipm++) {
          /* `c1' is the code of the permuted diagram `c' */
          for (c1 = 0, ipr = 0; ipr < npr; ipr++, ms++)
            if ((c >> ipr) & 1u) c1 |= *ms;
          if (m->map[c1] < 0) {
            m->map[c1] = (unqid_t) gid;
          } else die_if (m->map[c1] != gid,
            "%4d: error: corruption code: %#x %#x, graph %d, %d\n",
            inode, c, c1, m->map[c1], gid);
        }
        gid++;
      }
      free(masks);
      fprintf(stderr, "%4d: n %d, initialized, %d unique diagrams, %gs\n",
         inode, DG_N_, gid, 1.*(clock() - t0)/CLOCKS_PER_SEC);
      /* we set m->ng now, for everything is properly set now */
      m->ng = gid;
  END:
      ;
#ifdef _OPENMP
    } /* m->ng <= 0 */
  } /* omp critical */
#endif /* _OPENMP */
  return 0;
}



/* reset a dgmap_t, but don't free the point */
INLINE void dgmap_done(dgmap_t *m)
{
  if (m->ng > 0) {
    free(m->map);
    m->map = NULL;
    free(m->first);
    m->first = NULL;
    m->ng = 0;
  }
}



/* retrieve the diagram id */
INLINE unqid_t dg_getmapidx(const dg_t *g, dgword_t *c)
{
  DG_DEFN_(g)

  dgmap_init(&dgmap_[DG_N_], DG_N_);
  dg_encode(g, c);
  return dgmap_[DG_N_].map[*c];
}



/*  retrieve the diagram id */
INLINE unqid_t dg_getmapid(const dg_t *g)
{
  dgword_t c;
  return dg_getmapidx(g, &c);
}



/* the macro is slower than the direct version (due to dg_getmapid()) */
#define dg_biconnected_lookup(g) dg_biconnected_lookuplow(DG_GN_(g), dg_getmapid(g))

/* biconnectivity of unique diagrams
 * unique diagrams are few enough to allow private memories */
static int *dgmap_bc_[DGMAP_NMAX + 1] = {NULL};
#pragma omp threadprivate(dgmap_bc_)


/* check if a diagram is biconnected (lookup version)
 * use it only if `id' is known, otherwise it is slower */
INLINE int dg_biconnected_lookuplow(int n, unqid_t id)
{
  if (dgmap_bc_[DG_N_] == NULL) {
    /* initialize the look-up table */
    dg_t *g1 = dg_open(DG_N_);
    dgmap_t *m = dgmap_ + DG_N_;
    int k;

    dgmap_init(m, DG_N_);
    xnew(dgmap_bc_[DG_N_], m->ng);
    for (k = 0; k < m->ng; k++) {
      dgword_t c = m->first[k];
      dg_decode(g1, &c);
      dgmap_bc_[DG_N_][k] = dg_biconnected(g1);
    }
    dg_close(g1);
  }
  return dgmap_bc_[ DG_N_ ][ id ];
}
#endif /* defined(DGMAP_EXISTS) */



/* free all stock pointers */
INLINE void dgmap_free(void)
{
#ifdef DGMAP_EXISTS
  int k;

#pragma omp critical
  {
    for (k = 0; k <= DGMAP_NMAX; k++)
      dgmap_done(&dgmap_[k]);
  }

  /* dgmap_bc_[] is private */
  for (k = 0; k <= DGMAP_NMAX; k++)
    if (dgmap_bc_[k] != NULL)
      free(dgmap_bc_[k]);
#endif /* defined(DGMAP_EXISTS) */
}


#endif /* DGMAP_H__ */

