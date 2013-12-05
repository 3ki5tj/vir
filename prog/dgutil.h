#ifndef DGUTIL_H__
#define DGUTIL_H__



/* nonessential but useful routines */



/* choose a random edge, return the number of edges */
INLINE int dg_randedge(const dg_t *g, int *i0, int *i1)
{
  dgvs_t c, maski;
  dgword_t b;
  DGVS_DEFIQ_(iq)
  int i, j, k, ne, ipr, rr;
  DG_DEFN_(g)
  static int cnt[DG_NMAX];
#pragma omp threadprivate(cnt)

  *i0 = *i1 = 0; /* in case no edge exists */
  DGVS_CLEAR(maski);
  DGVS_ADD(maski, 0);
  for (ne = 0, i = 1; i < DG_N_; i++) {
    /* `maski' collects vertices with indices low than `i' */
    ne += cnt[i] = dg_degvs(g, i, maski);
    DGVS_ADD(maski, i) /* update `maski' */
  }
  rr = (int) (rnd0() * 2 * ne);
  ipr = rr / 2;
  /* go through pairs such that 0 <= j < i < N */
  DGVS_CLEAR(maski)
  DGVS_ADD(maski, 0)
  for (i = 1; i < DG_N_; ipr -= cnt[i], i++) {
    if (ipr < cnt[i]) { /* found it */
      DGVS_AND2(c, g->c[i], maski)
      for (k = -1, j = 0; j <= ipr; j++) {
        DGVS_FIRSTLOW(k, c, b, iq)
        DGVS_XOR1(c, b, iq) /* remove `b' from `c' */
      }
      die_if (k < 0, "i %d, k %d\n", i, k);
      /* (*i0, *i1) or (*i1, *i0) */
      if (rr % 2) { *i0 = i, *i1 = k; }
      else { *i0 = k, *i1 = i; }
      break;
    }
    DGVS_ADD(maski, i) /* update `maski' */
  }
  return ne;
}



/* construct `sg' by removing vertex `i0' from `g'
 * `sg' and `g' can be the same */
INLINE dg_t *dg_remove1(dg_t *sg, dg_t *g, int i0)
{
  int i, is = 0, n = g->n;
  dgvs_t maskl, maskh;
#if !DGVS_ONEWORD
  dgvs_t vs1, vs2;
#endif

#ifdef N
  die_if(n - 1 != N, "dg_remove1 cannot be used n %d with fixed N\n", n);
#endif
  DGVS_MKBITSMASK(maskl, i0)
  DGVS_NOT(maskh, maskl) /* maskh = ~maskl; */
  for (i = 0; i < n; i++) {
    if (i0 == i) continue;
#if DGVS_ONEWORD
    sg->c[is++] = (g->c[i] & maskl) | ((g->c[i] >> 1) & maskh);
#else /* !DGVS_ONEWORD */
    DGVS_AND2(vs1, g->c[i], maskl)
    DGVS_RSHIFT1(vs2, g->c[i])
    DGVS_AND(vs2, maskh)
    DGVS_OR2(sg->c[is], vs1, vs2)
    is++;
#endif /* DGVS_ONEWORD */
  }
  sg->n = n - 1;
  return sg;
}



#endif /* DGUTIL_H__ */

