#ifndef DGUTIL_H__
#define DGUTIL_H__

/* nonessential but useful routines */



/* choose a random edge, return the number of edges */
INLINE int dg_randedge(const dg_t *g, int *i0, int *i1)
{
  int i, j, k, ne, ipr, rr;
  DG_DEFN_(g)
  dgword_t c, b;
  static int cnt[DG_NMAX];
#pragma omp threadprivate(cnt)

  *i0 = *i1 = 0; /* in case no edge exists */
  for (ne = 0, i = 1; i < DG_N_; i++)
    ne += cnt[i] = bitcount(g->c[i] & MKBITSMASK(i));
  rr = (int) (rnd0() * 2 * ne);
  ipr = rr / 2;
  /* we will go through pairs such 0 <= j < i < N */
  for (i = 1; i < DG_N_; ipr -= cnt[i], i++) {
    if (ipr < cnt[i]) { /* found it */
      c = g->c[i] & MKBITSMASK(i);
      for (k = -1, j = 0; j <= ipr; j++, c ^= b) {
        BITFIRSTLOW(k, c, b);
      }
      die_if (k < 0, "i %d, k %d\n", i, k);
      if (rr % 2) { *i0 = i, *i1 = k; }
      else { *i0 = k, *i1 = i; }
      break;
    }
  }
  return ne;
}



/* construct `sg' by removing vertex `i0' from `g'
 * `sg' and `g' can be the same */
INLINE dg_t *dg_remove1(dg_t *sg, dg_t *g, int i0)
{
  int i, is = 0, n = g->n;
  dgword_t maskl, maskh;

#ifdef N
  die_if(n - 1 != N, "dg_remove1 cannot be used n %d with fixed N\n", n);
#endif
  maskl = MKBITSMASK(i0);
  maskh = ~maskl;
  for (i = 0; i < n; i++) {
    if (i0 == i) continue;
    sg->c[is++] = (g->c[i] & maskl) | ((g->c[i] >> 1) & maskh);
  }
  sg->n = n - 1;
  return sg;
}



#endif /* DGUTIL_H__ */

