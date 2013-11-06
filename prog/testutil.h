#ifndef TESTUTIL_H__
#define TESTUTIL_H__
/* utilities for test programs */
#include "dg.h"



#define dg_rndswitchedge(g, ned, nedmax, rnp) \
  dg_rndswitchedge0(g, ned, 0, 1, nedmax, rnp)

INLINE int dg_rndswitchedge0(dg_t *g, int *ned,
    int nedmin, double rnm, int nedmax, double rnp)
{
  int i = 0, j = 0, n = g->n, ipr, npr, acc;

  npr = n * (n - 1) / 2;
  /* randomly switch an edge */
  ipr = (int) (npr * rnd0());
  parsepairindex(ipr, n, &i, &j);
  if (dg_linked(g, i, j)) { /* unlink (i, j) */
    dg_unlink(g, i, j);
    acc = dg_biconnected(g);
    if ( acc && *ned <= nedmin )
      acc = (rnd0() < rnm);
    if ( acc ) {
      (*ned)--;
    } else { /* link back the edge */
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



/* adjust the rate of adding edges */
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


/* adjust the rate of removing edges */
INLINE void adjustrnm(int ned, int nedmin, int t, int nstadj,
    int *good, int *tot, double *rnm)
{
  *good += (ned >= nedmin);
  *tot += 1;
  if (t % nstadj == 0 && 1. * (*good) / (*tot) < 0.5) { /* adjust rnp */
    *rnm *= 0.5;
    printf("adjusting rnm to %g, ned %d, good %g%%\n",
        *rnm, ned, 100.*(*good)/(*tot));
    *good = *tot = 0;
  }
}
#endif
