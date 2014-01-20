#ifndef DGRJW_H__
#define DGRJW_H__
/* compute the overall weight of a configuration
 * by the method of Richard J. Wheatley, PRL 110, 200601 (2013)
 * iterative, small memory version 1 */
#include "dgcsep.h"
#include "dgsc.h"
#include <time.h>
#include <limits.h>



/* Since RJWNMAX is currently much less than DG_WORDBITS,
 * we assume a one-word vertex set, even if DGVS_ONEWORD == 0 */

#if !defined(RJW64) && !defined(RJWDBL) && \
    (defined(RJW32) || (defined(N) && N <= 14 && DG_WORDBITS == 32))
  typedef int32_t dgrjw_fb_t;
  /* for n <= 14, max |fb(n)| = 12! = 479001600 < |DGRJW_FBDIRTY|
   * Note the code *may* fail before this limit, use testrjw.c to test
   * Using the 32-bit integer is advantageous for saves memory
   * and runs faster on a 32-bit machine */
  #ifndef RJWNMAX
  #define RJWNMAX 14
  #endif

#elif !defined(RJWDBL) && (!defined(N) || N <= 22) /* 64-bit RJW */
  /* this is the default branch for with N > 22
   * there is usually a memory problem */

  typedef int64_t dgrjw_fb_t;
  /* for n <= 22, max |fb(n)| = 20! = 2.4e18 < |DGRJW_FBDIRTY| */
  /* Note: the threshold is based on the limit 20! = 2.4e18
   *  and 21! cannot be contained in a 64-bit integer
   * The code may however fail before that due to the cancellation
   *  in intermediate steps
   * Another limit is that memory > 2^(n + 1) */
  #ifndef RJWNMAX
  #define RJWNMAX 22
  #endif /* !defined(RJWNMAX) */

#else /* double RJW */

  #ifndef DGRJW_DOUBLE
  #define DGRJW_DOUBLE
  #endif /* !defined(DGRJW_DOUBLE) */

  typedef double dgrjw_fb_t;

  #ifndef RJWNMAX
  #define RJWNMAX 26
  #endif

#endif /* RJW32 or RJW64 */



/* articulation point
 * an unsigned char is sufficient for N < 256 */
typedef unsigned char dgrjw_ap_t;



/* obtain a list of numbers 0..2^n - 1 by the number of bits */
INLINE void dgrjw_idbybits(int n, unsigned *arr)
{
  int m, top, v, id = 0, st[DG_NMAX + 1];
  unsigned x; /* the number represented by the stack */

  arr[id++] = 0;
  for (m = 1; m <= n; m++) {
    x = 0;
    /* enumerate numbers from 1 to 2^n with m bits */
    st[top = 0] = 0;
    while (top >= 0) {
      //printf("n %d, m %d, top %d, st[top] %d\n", n, m, top, st[top]);
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



/* compute the Boltzmann weight, the sum of all diagrams
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dgrjw_fq(dgvs_t *c, dgword_t vs)
{
  dgword_t w, b;

  /* if there is a bond, i.e., r(i, j) < 1, then Eq = 0 */
  for (w = vs; w; w ^= b) {
    /* if c[i] share vertices with vs, there is bond
     * the Boltzmann weight = \prod_(ij) e_ij, and it allows no clash
     * therefore return zero immediately */
    if (DGVS_FIRSTWORD(c[bitfirstlow(w, &b)]) & vs) return 0;
  }
  return 1;
}



/* increment ms1 in the limits of bits of (ms1 | ms2) */
#define DGRJW_INC(ms1, ms2) DGRJW_INCa(ms1, ms2)

/* version 1 */
#define DGRJW_INCa(ms1, ms2) { dgword_t nms2_, lbit_, hbits_; \
  nms2_ = (dgword_t) (-(ms2)); /* -ms2 and ms2 share bits higher than the lowest bit of ms2 */ \
  lbit_ = (ms2) & nms2_; /* the & yields the lowest bit lowbit_ */ \
  hbits_ = (ms2) ^ nms2_; /* the collections of bits higher than lbit_ */ \
  (ms1) &= hbits_; /* this wipes out bits lower or equal to lbit_ */ \
  (ms1) |= lbit_;  /* this add the lowest bit */ }

/* version 2 */
#define DGRJW_INCb(ms1, ms2) { dgword_t lbit_; \
    lbit_ = (ms2) & (-ms2); /* find the lowest empty bit (first unused variable vertex) */ \
    (ms1) ^= lbit_; /* add this bit */ \
    (ms1) &= ~(lbit_ - 1); /* clear all lower bits */ }

/* version 3 */
#define DGRJW_INCc(ms1, ms2) { int i = BITFIRSTNZ(ms2); \
  ms1 >>= i; /* clear the lower bits by shifting i bits */ \
  ms1 ^= 1; /* add the lowest bit, `^' can be replaced by `+' or `|' with little difference */ \
  ms1 <<= i; /* shift back */ }



/* compute the Boltzmann weight (fq) and the sum of connected diagrams (fc)
 * also a few values for fa and fb
 * fqarr extends to faarr, fcarr extends to fbarr
 * by Wheatley's recursion formula, for all diagrams
 * `c' is the connectivity matrix */
INLINE void dgrjw_prepare(const dg_t *g,
    dgrjw_fb_t *fcarr, dgrjw_fb_t *fqarr,
    unsigned *idbybits, dgrjw_ap_t *aparr)
{
  dgrjw_fb_t fc;
  dgword_t vs, vs1, ms, ms1, ms2, b1, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */

  /* we loop from smaller subsets to larger ones */
  /* diagrams with zero or one vertex */
  for (id = 0; id <= (unsigned) DG_N_; id++) {
    vs = idbybits[id];
    die_if (bitcount(vs) > 1, "n %d, id %d, bad count(%#x) = %d\n", DG_N_, id, vs, bitcount(vs));
    fcarr[vs] = fqarr[vs] = 1;
  }

  /* diagrams with two vertices */
  for (; id <= (unsigned) DG_N_ * (DG_N_ + 1) / 2; id++) {
    int fq;
    vs = idbybits[id];
    die_if (bitcount(vs) != 2, "n %d, id %d, bad count(%#x) = %d\n", DG_N_, id, vs, bitcount(vs));
    fq = dgrjw_fq(g->c, vs);
    fc = -!fq;
    fqarr[vs] = fq;
    fcarr[vs] = fc;
  }
  /* everything is known for diagrams of two or fewer vertices */

  /* diagrams with three or more vertices */
  for (; id < vsmax; id++) {
    vs = idbybits[id];
    fc = dgrjw_fq(g->c, vs);
    fqarr[vs] = fc; /* start with the Boltzmann weight */
    if ( !dg_connectedvs(g, vs) ) { /* disconnected subset */
      fcarr[vs] = 0;
      aparr[vs] = 0; /* fb(vs) == 0 */
    } else { /* connected diagram */
      b1 = vs & (-vs);
      ms = vs ^ b1;
      /* loop over proper subsets of `vs'
       * the first component contains `b1' */
      for (ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
        fc -= fcarr[ms1 ^ b1] * fqarr[ms2];
        DGRJW_INCa(ms1, ms2); /* update the subset `ms1' */
      }
      fcarr[vs] = fc;

      /* check articulation points */
      aparr[vs] = (dgrjw_ap_t) (DG_N_ + 1);
      for ( vs1 = vs; vs1 != 0; ) {
        b1 = vs1 & (-vs1);
        vs1 ^= b1;
        if ( !dg_connectedvs(g, vs ^ b1) ) {
          /* fb(vs, v) == 0 for v = BIT2ID(b1) + 1 to n - 1 */
          aparr[vs] = (dgrjw_ap_t) (BIT2ID(b1) + 1);
          break;
        }
      }
    } /* connected / disconnected diagrams */
  }
}



/* return fb, assuming dgrjw_prepare() has been called */
INLINE dgrjw_fb_t dgrjw_iter(const dg_t *g,
    dgrjw_fb_t *fbarr, unsigned *idbybits, dgrjw_ap_t *aparr)
{
  int v, iold, inew;
  size_t idold, idnew, jdold, jdnew;
  dgrjw_fb_t fa;
  dgword_t vs, vsnew, ms, ms1, ms2, b1, bv, b1v, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */

  /* one- and two-vertex subsets, biconnected == connected */
  for ( id = 0; id <= (dgword_t) (DG_N_ * (DG_N_ + 1) / 2); id++ ) {
    vs = idbybits[id];
    fbarr[ vsmax | vs ] = fbarr[ vs ];
  }

  /* loop over the position of the first clique separator */
  for ( v = 0; v < DG_N_; v++ ) {
    iold = v % 2;
    inew = (v + 1) % 2;
    idold = iold * vsmax;
    idnew = inew * vsmax;

    bv = MKBIT(v);

    /* three-vertex subsets */
    for ( id = DG_N_ * (DG_N_ + 1) / 2 + 1 ; id < vsmax; id++ ) {
      vs = idbybits[id];
      vsnew = idnew + vs;

      /* compute fa and fb */
      if ( aparr[vs] <= v ) {
        fbarr[vsnew] = 0;
      } else {
        fa = 0;
        /* skip if the vertex set does not contain v */
        if ( (bv & vs) != 0 ) {
          ms = vs ^ bv;
          b1 = ms & (-ms); /* lowest vertex */
          b1v = b1 ^ bv;
          ms ^= b1; /* remove the fixed vertices `b1' and `bv' from `vs' */

          jdnew = idnew | b1v;
          jdold = idold | bv;
          for ( ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
            fa += fbarr[jdnew | ms1] * fbarr[jdold | ms2];
            /* update the subset `ms1' */
            DGRJW_INCa(ms1, ms2);
          }
        }
        fbarr[vsnew] = fbarr[idold | vs] - fa; /* set fb */
      }
    }
  }
  return fbarr[vsmax * (DG_N_ % 2) + vsmax - 1];
}



static int dgrjw_nmax_;
static dgrjw_fb_t *dgrjw_fbarr_;
/* every thread needs its own memory to do independent calculation
 * so these variables must be thread private */
#pragma omp threadprivate(dgrjw_nmax_, dgrjw_fbarr_)

static unsigned *dgrjw_idbybits_;
static int dgrjw_idn_ = 0;
#pragma omp threadprivate(dgrjw_idbybits_, dgrjw_idn_)

static dgrjw_ap_t *dgrjw_aparr_;
#pragma omp threadprivate(dgrjw_aparr_)



/* compute the sum of biconnected diagrams by Wheatley's method
 * This is a low level function and the test of clique separator
 * is not done here. */
INLINE dgrjw_fb_t dgrjw_fb(const dg_t *g)
{
  DG_DEFN_(g)
  size_t size = 0;

  /* the memory requirement is 2^(n + 1) * (n + 1) * sizeof(dgrjw_fb_t) */
  if (dgrjw_fbarr_ == NULL) {
    /* even if we know N, we still use dynamic allocation
     * because this function might not be called at all
     * in high dimensions */
    dgrjw_nmax_ = DG_N_;
    size = (size_t) 1u << dgrjw_nmax_; /* 2^nmax */
    xnew(dgrjw_fbarr_, size * 2);
    xnew(dgrjw_idbybits_, size);
    xnew(dgrjw_aparr_, size);
    fprintf(stderr, "%4d: dgrjw allocated %gMB memory for n %d\n", inode,
        (2. * sizeof(dgrjw_fb_t) + sizeof(unsigned) + sizeof(dgrjw_ap_t))
        * size / (1024*1024), dgrjw_nmax_);
  }
#ifndef N /* if N is fixed, no need to reallocate */
  else if (DG_N_ > dgrjw_nmax_) {
    dgrjw_nmax_ = DG_N_;
    size = (size_t) 1u << dgrjw_nmax_; /* 2^nmax */
    xrenew(dgrjw_fbarr_, size * 2);
    xrenew(dgrjw_idbybits_, size);
    xrenew(dgrjw_aparr_, size);
    fprintf(stderr, "%4d: dgrjw reallocated %gMB memory for n %d\n", inode,
        (2. * sizeof(dgrjw_fb_t) + sizeof(unsigned) + sizeof(dgrjw_ap_t))
        * size / (1024*1024), dgrjw_nmax_);
  }
#endif

#if !defined(N) || N < 32
  size = ((size_t) 1u << DG_N_); /* 2^n */
#endif

  if (dgrjw_idn_ != DG_N_) { /* re-sort diagrams */
    dgrjw_idbybits(DG_N_, dgrjw_idbybits_);
    dgrjw_idn_ = DG_N_;
  }

  /* pre-compute all fc and fq values */
  dgrjw_prepare(g, dgrjw_fbarr_, dgrjw_fbarr_ + size,
                   dgrjw_idbybits_, dgrjw_aparr_);
  return dgrjw_iter(g, dgrjw_fbarr_, dgrjw_idbybits_, dgrjw_aparr_);
}



#define dg_fb(g) \
  dg_fb0(g, DGCSEP_DEFAULTMETHOD, NULL, NULL)

/* directly compute the sum of biconnected diagrams by various strategies
 * a mixed strategy of dgrjw_fb() and dg_sc()
 * nocsep = 1 means that the graph has been tested with no clique separator
 *        = 0 means it MAY have clique separators
 * *ned: number of edges; degs: degree sequence
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_fb0(const dg_t *g, int csepmethod, int *ned, int *degs)
{
  int err, nedges = -1;
  DG_DEFN_(g)
  double sc;

  if ( ned == NULL ) ned = &nedges;
  sc = dgsc_spec0(g, csepmethod, ned, degs, &err);
  if ( err == 0 ) {
    return DG_SC2FB(sc, *ned);
  } else if ( *ned <= 2*DG_N_ - 3 || DG_N_ > RJWNMAX) {
    return DG_SC2FB(dgsc_do(g), *ned);
  } else { /* hs_rjw() requires 2^(n + 1) * (n + 1) memory */
    return (double) dgrjw_fb(g);
  }
}



/* free all stock objects */
INLINE void dgrjw_free(void)
{
  if (dgrjw_fbarr_ != NULL) {
    free(dgrjw_fbarr_);
    dgrjw_fbarr_ = NULL;
  }
  if (dgrjw_idbybits_ != NULL) {
    free(dgrjw_idbybits_);
    dgrjw_idbybits_ = NULL;
  }
  if (dgrjw_aparr_ != NULL) {
    free(dgrjw_aparr_);
    dgrjw_aparr_ = NULL;
  }
  dgrjw_nmax_ = 0;
}



#endif

