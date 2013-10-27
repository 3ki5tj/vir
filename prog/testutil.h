#ifndef TESTUTIL_H__
#define TESTUTIL_H__
/* utilities for test programs */
#include "dg.h"


INLINE int dg_rndswitchedge(dg_t *g, int *ned, int nedmax, double rnp)
{
  int i, j, n = g->n, ipr, npr, acc;

  npr = n * (n - 1) / 2;
  /* randomly switch an edge */
  ipr = (int) (npr * rnd0());
  parsepairindex(ipr, n, &i, &j);
  if (dg_linked(g, i, j)) { /* unlink (i, j) */
    dg_unlink(g, i, j);
    acc = dg_biconnected(g);
    if ( acc ) {
      (*ned)--;
    } else {
      dg_link(g, i, j);
    }
  } else { /* link (i, j) */
    acc = 1;
    if (*ned >= nedmax) /* avoid increasing edges */
      acc = (rnd0() < rnp);
    if (acc) {
      (*ned)++;
      dg_link(g, i, j);
    }
  }
  return acc;
}



INLINE void dg_linkpairs(dg_t *g, int pair[][2])
{
  int i;

  for (i = 0; pair[i][0] >= 0; i++)
    dg_link(g, pair[i][0], pair[i][1]);
}



/* adjust the rate of increasing edges */
INLINE void adjustrnp(int ned, int nedmax, int t, int nstadj,
    int *good, int *tot, double *rnp)
{
  *good += (ned <= nedmax);
  *tot += 1;
  if (t % nstadj == 0 && 1. * (*good) / (*tot) < 0.5) { /* adjust rnp */
    *rnp *= 0.5;
    printf("adjusting rnp to %g, ned %d, good %g%%\n",
        *rnp, ned, 100.*(*good)/(*tot));
    *good = *tot = 0;
  }
}
#endif
