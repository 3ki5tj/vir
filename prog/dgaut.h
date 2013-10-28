#ifndef DGAUT_H__
#define DGAUT_H__
/* canonical label */

#include "dg.h"


/* set nauty parameteters */
#define WORDSIZE        CODEBITS
#define MAXN            DG_NMAX
#define ONE_WORD_SETS
#include "nau0s.h"



/* canonical label */
INLINE dg_t *dg_canlabel(dg_t *gout, const dg_t *gin)
{
  static code_t g0[DG_NMAX], gc[DG_NMAX];
  static int lab[DG_NMAX], ptn[DG_NMAX], orbits[DG_NMAX];
  static DEFAULTOPTIONS_GRAPH(options);
  static statsblk stats;
#pragma omp threadprivate(g0, gc, lab, ptn, orbits, options, stats)
  int i, n = gin->n;

  /* nauty uses the highest bit for the first index */
  for (i = 0; i < n; i++) g0[i] = bitreverse(gin->c[i]);
  options.getcanon = TRUE;
  densenauty(g0, lab, ptn, orbits, &options, &stats, 1, n, gc);

  gout->n = n;
  for (i = 0; i < n; i++) gout->c[i] = bitreverse(gc[i]);
  return gout;
}



#endif /* DGAUT_H__ */

