#ifndef DGUTIL_H__
#define DGUTIL_H__



/* nonessential but useful routines */



#include "dg.h"



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



/* table of factorials */
static double dg_facts_[DG_NMAX + 1] = {0};
static int dg_facts_init_ = 0;
#pragma omp threadprivate(dg_facts_, dg_facts_init_)

#define DG_CALC_FACTS() if ( !dg_facts_init_ ) { \
  dg_facts_[0] = dg_facts_init_ = 1; \
  for (i = 1; i <= DG_NMAX; i++) \
    dg_facts_[i] = dg_facts_[i - 1] * i; \
}



/* table of 2^n */
static double dg_pow2s_[DG_NMAX * DG_NMAX + 1] = {0};
static int dg_pow2s_init_ = 0;
#pragma omp threadprivate(dg_pow2s_, dg_pow2s_init_)

#define DG_CALC_POW2S() if ( !dg_pow2s_init_ ) { \
  dg_pow2s_[0] = dg_pow2s_init_ = 1; \
  for (i = 1; i <= DG_NMAX * DG_NMAX; i++) \
    dg_pow2s_[i] = dg_pow2s_[i - 1] * 2; \
}



/* obtain a list of numbers 0..2^n - 1 by the number of bits */
INLINE void dg_idbybits(int n, unsigned *arr)
{
  int m, top, v, id = 0, st[DG_NMAX + 1];
  unsigned x; /* the number represented by the stack */

  arr[id++] = 0;
  for (m = 1; m <= n; m++) {
    x = 0;
    /* enumerate numbers from 1 to 2^n with m bits */
    st[top = 0] = 0;
    while (top >= 0) {
      if ( st[top] < n ) { /* push */
        v = st[top++];
        x ^= 1u << v;
        /* push if we still have numbers to fill */
        if (top < m) {
          /* the number on the next level must be greater than
           * the number on this level */
          st[top] = v + 1;
          continue;
        } else { /* we are at the top, print and fall through to pop */
          arr[id++] = x;
        }
      }
      /* pop */
      if (--top >= 0) {
        v = st[top];
        x ^= 1u << v; /* st[top] is no longer available on top */
        st[top] = v + 1;
      }
    }
  } /* end of the round for m bits */
}



/* prepare an idbybits array for n,
 * `*n0' is the current n, `*nmax' is the capacity */
INLINE unsigned *dg_prep_idbybits(int n, unsigned *arr,
    int *n0, int *nmax)
{
  size_t size = (size_t) 1 << n;

  if ( arr == NULL ) {
    xnew(arr, size);
    *nmax = n; *n0 = -1;
  } else if ( n > *nmax ) {
    xrenew(arr, size);
    *nmax = n; *n0 = -1;
  }
  if ( n != *n0 ) {
    dg_idbybits(n, arr);
    *n0 = n;
  }
  return arr;
}



/* increment ms1 in the limits of bits of (ms1 | ms2) */
#define DGVS_INC(ms1, ms2) DGVS_INCa(ms1, ms2)

/* version 1 */
#define DGVS_INCa(ms1, ms2) { dgword_t nms2_, lbit_, hbits_; \
  nms2_ = (dgword_t) (-(ms2)); /* -ms2 and ms2 share bits higher than the lowest bit of ms2 */ \
  lbit_ = (ms2) & nms2_; /* the & yields the lowest bit lowbit_ */ \
  hbits_ = (ms2) ^ nms2_; /* the collections of bits higher than lbit_ */ \
  (ms1) &= hbits_; /* this wipes out bits lower or equal to lbit_ */ \
  (ms1) |= lbit_;  /* this add the lowest bit */ }

/* version 2 */
#define DGVS_INCb(ms1, ms2) { dgword_t lbit_; \
    lbit_ = (ms2) & (-ms2); /* find the lowest empty bit (first unused variable vertex) */ \
    (ms1) ^= lbit_; /* add this bit */ \
    (ms1) &= ~(lbit_ - 1); /* clear all lower bits */ }

/* version 3 */
#define DGVS_INCc(ms1, ms2) { int i = BITFIRSTNZ(ms2); \
  ms1 >>= i; /* clear the lower bits by shifting i bits */ \
  ms1 ^= 1; /* add the lowest bit, `^' can be replaced by `+' or `|' with little difference */ \
  ms1 <<= i; /* shift back */ }



#endif /* DGUTIL_H__ */

