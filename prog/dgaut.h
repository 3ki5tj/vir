#ifndef DGAUT_H__
#define DGAUT_H__
/* canonical label */

#include "dg.h"


/* set parameteters for the program `nauty' (No AUTomorphism, Yes?) */
#define WORDSIZE        CODEBITS
#define MAXN            DG_NMAX
#define ONE_WORD_SETS 1 /* use one-word set when possible */
#include "nau0s.h"



/* canonical label */
INLINE dg_t *dg_canlabel(dg_t *gout, const dg_t *gin)
{
  static code_t g0[DG_NMAX], gc[DG_NMAX];
  static int lab[DG_NMAX], ptn[DG_NMAX], orbits[DG_NMAX];
  static DEFAULTOPTIONS_GRAPH(options);
  static statsblk stats;
#pragma omp threadprivate(g0, gc, lab, ptn, orbits, options, stats)
  int i;
  DG_DEFN_(gin);

  /* nauty uses the highest bit for the first index */
  for (i = 0; i < DG_N_; i++) g0[i] = bitreverse(gin->c[i]);
  options.getcanon = TRUE;
  densenauty(g0, lab, ptn, orbits, &options, &stats, 1, DG_N_, gc);

  gout->n = DG_N_;
  for (i = 0; i < DG_N_; i++) gout->c[i] = bitreverse(gc[i]);
  return gout;
}



#endif /* DGAUT_H__ */

