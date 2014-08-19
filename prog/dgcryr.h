#ifndef DGCRYR_H__
#define DGCRYR_H__
/* compute the direct correlation function
 * and the cavity function */
#include "dg.h"



/* iteratively compute the signed star content for y(r)
 * which is computed as follows
 * First, join the root points a-b, if they are not adjacent in g
 * Then, the signed star content of a-b is equal to
 * the number of biconnected subgraphs of g
 * with even number of missing edges (except a-b)
 * less the number of biconnected subgraphs of g
 * with odd number of missing edges (except a-b)
 * On return g should be the same */
__inline double dgsc_yiter(dg_t *g, double *nr, int a, int b)
{
  static int ed[DG_NMAX * (DG_NMAX - 1)/2][2]; /* edges */
  static int st[DG_NMAX * (DG_NMAX - 1)/2 + 2]; /* state */
#pragma omp threadprivate(ed, st)
  int n = g->n, top, ied, med, vi, vj, degi, degj, eab;
  double fb, ring;
  int ned = -1, degs[DG_NMAX];

  if ( a > b ) vi = a, a = b, b = vi;

  die_if (b < 0, "(%d-%d) is not a valid edge\n", a, b);

  /* connect a-b */
  if ( !dg_linked(g, a, b) ) {
    dg_link(g, a, b);
    eab = 1;
  } else {
    eab = 0;
  }

  /* collect edges */
  ned = dg_degs(g, degs);
  for (med = 0, vi = 0; vi < n - 1; vi++)
    if ( degs[vi] > 2 )
      for (vj = vi + 1; vj < n; vj++)
        if ( degs[vj] > 2 && dg_linked(g, vi, vj) )
          ed[med][0] = vi, ed[med][1] = vj, med++;
  if (nr == NULL) nr = &ring;

  st[ top = 0 ] = -1;
  fb = 1 - ned % 2 * 2;
  *nr = (ned == n) ? 1 : 0;
  while (top >= 0) {
    /* search over remaining edges */
    for (ied = st[top] + 1; ied < med; ied++) {
      /* removing an edge adjacent to a degree-2 vertex
       * makes a graph not biconnected */
      if ( (degi = degs[ vi = ed[ied][0] ]) < 3
        || (degj = degs[ vj = ed[ied][1] ]) < 3 
        || (vi == a && vj == b) )
        continue;
      dg_unlink(g, vi, vj);
      if ( dg_biconnected(g) ) { /* push */
        degs[vi] = degi - 1, degs[vj] = degj - 1;
        st[top] = ied;
        --ned;
        fb += 1 - ned % 2 * 2;
        if (ned == n) *nr += 1;
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
      ned++;
    }
  }
  
  if ( eab ) dg_unlink(g, a, b); /* disconnect a-b */
  fb *= -1; /* remove the bond between a and b */
  
  return fb;
}





/* histogram for correlation functions */
typedef struct {
  double xmin;
  double xmax;
  double dx;
  double sphr;
  double vol;
  double *arr;
  int dim;
  int order;
  int npt; /* number of points along r */
  unsigned long nsampi;
  double nsampd;
} hscr_t;



INLINE hscr_t *hscr_open(real xmax, real dx, int dim, int order)
{
  hscr_t *hs;
  int i;

  xnew(hs, 1);
  hs->xmin = 0;
  hs->xmax = xmax;
  hs->dx = dx;
  hs->dim = dim;
  hs->order = order;
  hs->npt = (int)(xmax/dx + 0.99999999);
  hs->nsampi = 0;
  hs->nsampd = 0;
  hs->sphr = (dim % 2 + 1) * dim;
  for (i = 2 + dim % 2; i <= dim; i += 2 )
    hs->sphr *= M_PI * 2 / i;
  hs->vol = hs->sphr / dim;
  xnew(hs->arr, hs->npt);
  return hs;
}



/* increase the count of the histogram */
INLINE void hscr_inc(hscr_t *hs, int cnt)
{
  if ( (hs->nsampi += cnt) > 1000000000 ) {
    hs->nsampd += (double) hs->nsampi;
    hs->nsampi = 0;
  }
}



/* add an entry, but do not update the sample size */
INLINE void hscr_add0(hscr_t *hs, double x, double fb)
{
  int ix = (int)((x - hs->xmin)/hs->dx);
  if (ix < hs->npt) hs->arr[ix] += fb;
}



/* add an entry and update the sample size */
INLINE void hscr_add(hscr_t *hs, double x, double fb)
{
  hscr_add0(hs, x, fb);
  hscr_inc(hs, 1);
}



INLINE int hscr_save(hscr_t *hs, const char *fn, double norm)
{
  FILE *fp;
  int i, imax;
  double y, r;

  /* Zn/V^(n-1) --> Zn */
  norm *= pow(hs->vol, hs->order - 1);
  /* divide it by (n - 2)! for permutations of the black dots */
  for ( i = 2; i <= hs->order - 2; i++ ) norm /= i;
  for ( imax = hs->npt - 1; imax >= 0; imax-- )
    if ( fabs(hs->arr[imax]) > 0 ) break;
  if ( imax < 0 ) return -1;
  /* flush the counter */
  hs->nsampd += (double) hs->nsampi;
  hs->nsampi = 0;
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d %d %.0f %g %g %g %g\n",
      hs->dim, hs->order, imax + 1,
      hs->nsampd, hs->dx, norm, hs->sphr, hs->vol);
  for ( i = 0; i <= imax; i++ ) {
    r = (i + .5) * hs->dx;
    y = hs->arr[i]/(hs->nsampd * hs->sphr * pow(r, hs->dim - 1) * hs->dx);
    if (norm > 0) y *= norm;
    fprintf(fp, "%12.7f %22.14e\n", r, y);
  }
  fclose(fp);
  return 0;
}



INLINE void hscr_close(hscr_t *hs)
{
  if (!hs) return;
  if (hs->arr) free(hs->arr);
  memset(hs, 0, sizeof(*hs));
}



#endif /* defined DGCRYR_H__ */

