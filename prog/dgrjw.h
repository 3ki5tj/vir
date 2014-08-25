#ifndef DGRJW_H__
#define DGRJW_H__
/* compute the overall weight of a configuration
 * by the method of Richard J. Wheatley, PRL 110, 200601 (2013)
 * iterative, small memory version 1 */
#include "dgutil.h"
#include "dgcsep.h"
#include "dgsc.h"
#include "dgring.h"
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
  /* this is the default branch for undefined N and N <= 22 */

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



/* compute the Boltzmann weight, the sum of all diagrams
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dgrjw_fq(dgvs_t *c, dgword_t vs)
{
  dgword_t w, b;

  /* if there is a bond, i.e., r(i, j) < 1, then fq = 0 */
  for (w = vs; w; w ^= b) {
    /* if c[i] share vertices with vs, there is bond
     * the Boltzmann weight = \prod_(ij) e_ij, and it allows no clash
     * therefore return zero immediately */
    if (DGVS_FIRSTWORD(c[bitfirstlow(w, &b)]) & vs) return 0;
  }
  return 1; /* (vs != 0); */
}



/* compute the Boltzmann weight (fq) and the sum of connected diagrams (fc)
 * also a few values for fa and fb
 * `fqarr' extends to `faarr', `fcarr' extends to `fbarr'
 * `vsbysize[i]' gives the ith subset of `g->n' vertices,
 *    where the subsets are sorted by the number of vertices contained
 * `aparr[vs]' is the index array of the lowest articulation point + 1
 *    in the diagram induced by the vertex set `vs'
 * by Wheatley's recursion formula, for all diagrams */
static void dgrjw_prepare(const dg_t *g,
    dgrjw_fb_t *fcarr, dgrjw_fb_t *fqarr,
    unsigned *vsbysize, dgrjw_ap_t *aparr)
{
  dgrjw_fb_t fc;
  dgword_t vs, vs1, ms, ms1, ms2, b1, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */

  /* we loop from smaller subsets to larger ones */
  /* diagrams with zero or one vertex */
  for (id = 0; id <= (unsigned) DG_N_; id++) {
    vs = vsbysize[id];
    die_if (bitcount(vs) > 1, "n %d, id %d, bad count(%#x) = %d\n", DG_N_, id, vs, bitcount(vs));
    fcarr[vs] = fqarr[vs] = 1;
  }

  /* diagrams with two vertices */
  for (; id <= (unsigned) DG_N_ * (DG_N_ + 1) / 2; id++) {
    vs = vsbysize[id];
    die_if (bitcount(vs) != 2, "n %d, id %d, bad count(%#x) = %d\n", DG_N_, id, vs, bitcount(vs));
    fqarr[vs] = dgrjw_fq(g->c, vs);
    fcarr[vs] = -!fqarr[vs];
  }

  /* diagrams with three or more vertices */
  for (; id < vsmax; id++) {
    vs = vsbysize[id];

    fc = dgrjw_fq(g->c, vs);
    fqarr[vs] = fc; /* start with the Boltzmann weight */

    if ( !dg_connectedvs(g, DGVS_W2VS(vs)) ) { /* disconnected subset */
      fcarr[vs] = 0;
      aparr[vs] = 0; /* fb(vs) == 0 */
    } else { /* connected diagram */
      b1 = vs & (-vs);
      ms = vs ^ b1;
      /* loop over proper subsets of `vs'
       * the first component contains `b1' */
      for (ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
        fc -= fcarr[ms1 ^ b1] * fqarr[ms2];
        DGVS_INCa(ms1, ms2); /* update the subset `ms1' */
      }
      fcarr[vs] = fc;

      /* check articulation points */
      aparr[vs] = (dgrjw_ap_t) (DG_N_ + 1);
      for ( vs1 = vs; vs1 != 0; ) {
        b1 = vs1 & (-vs1);
        vs1 ^= b1;
        if ( !fcarr[vs ^ b1] ) { // if ( !dg_connectedvs(g, vs ^ b1) ) {
          /* fb(vs, v) == 0 for v = BIT2ID(b1) + 1 to n - 1 */
          aparr[vs] = (dgrjw_ap_t) (BIT2ID(b1) + 1);
          break;
        }
      }
    } /* connected / disconnected diagrams */
  }
}



/* return fb, assuming dgrjw_prepare() has been called */
static dgrjw_fb_t dgrjw_iter(const dg_t *g,
    dgrjw_fb_t *fbarr, unsigned *vsbysize, dgrjw_ap_t *aparr)
{
  int v, iold, inew;
  size_t idold, idnew, jdold, jdnew;
  dgrjw_fb_t fa;
  dgword_t vs, vsnew, ms, ms1, ms2, b1, bv, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */

  /* one- and two-vertex subsets, biconnected == connected */
  for ( id = 0; id <= (dgword_t) (DG_N_ * (DG_N_ + 1) / 2); id++ ) {
    vs = vsbysize[id];
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
      vs = vsbysize[id];
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
          ms ^= b1; /* remove the fixed vertices `b1' and `bv' from `vs' */

          jdnew = idnew | b1 | bv;
          jdold = idold | bv;
          for ( ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
            fa += fbarr[jdnew | ms1] * fbarr[jdold | ms2];
            /* update the subset `ms1' */
            DGVS_INCa(ms1, ms2);
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

/* subsets of n vertices sorted by the size */
static unsigned *dgrjw_vsbysize_;
static int dgrjw_vsn_ = -1, dgrjw_vsnmax_ = -1;
#pragma omp threadprivate(dgrjw_vsbysize_, dgrjw_vsn_, dgrjw_vsnmax_)

static dgrjw_ap_t *dgrjw_aparr_;
#pragma omp threadprivate(dgrjw_aparr_)



/* compute the sum of biconnected diagrams by Wheatley's method
 * This is a low level function and the test of clique separator
 * is not done here.
 * The memory requirement is 2^(n-20)*21MB */
static dgrjw_fb_t dgrjw_fb(const dg_t *g)
{
  DG_DEFN_(g)
  size_t size = 0;

#if !defined(N) || N < 32
  size = ((size_t) 1u << DG_N_); /* 2^n */
#endif

  /* the memory requirement is 2^(n + 1) * (n + 1) * sizeof(dgrjw_fb_t) */
  if (dgrjw_fbarr_ == NULL) {
    /* even if we know N, we still use dynamic allocation
     * because this function might not be called at all
     * in high dimensions */
    dgrjw_nmax_ = DG_N_;
    xnew(dgrjw_fbarr_, size * 2);
    xnew(dgrjw_aparr_, size);
    fprintf(stderr, "%4d: dgrjw allocated %gMB memory for n %d\n", inode,
        (2. * sizeof(dgrjw_fb_t) + sizeof(unsigned) + sizeof(dgrjw_ap_t))
        * size / (1024*1024), dgrjw_nmax_);
  } else if (DG_N_ > dgrjw_nmax_) {
    dgrjw_nmax_ = DG_N_;
    xrenew(dgrjw_fbarr_, size * 2);
    xrenew(dgrjw_aparr_, size);
    fprintf(stderr, "%4d: dgrjw reallocated %gMB memory for n %d\n", inode,
        (2. * sizeof(dgrjw_fb_t) + sizeof(unsigned) + sizeof(dgrjw_ap_t))
        * size / (1024*1024), dgrjw_nmax_);
  }

  dgrjw_vsbysize_ = dg_prep_vsbysize(DG_N_, dgrjw_vsbysize_,
      &dgrjw_vsn_, &dgrjw_vsnmax_);

  /* pre-compute all fc and fq values */
  dgrjw_prepare(g, dgrjw_fbarr_, dgrjw_fbarr_ + size,
                   dgrjw_vsbysize_, dgrjw_aparr_);
  return dgrjw_iter(g, dgrjw_fbarr_, dgrjw_vsbysize_, dgrjw_aparr_);
}



#define dg_fb(g) dg_fbnr(g, NULL)
#define dg_fb0(g, method, csepmethod, ned, degs) \
    dg_fbnr0(g, NULL, method, csepmethod, ned, degs)
#define dg_fbnr(g, nr) dg_fbnr0(g, nr, \
    DGSC_DEFAULTMETHOD, DGCSEP_DEFAULTMETHOD, NULL, NULL)

/* compute the sum of biconnected diagrams by various strategies
 * a mixed strategy of dgrjw_fb() and dgsc_fb()
 * `method': method of directly computing the star content
 * `csepmethod': method of detecting clique separators
 *    0 means not to detect clique separators
 * *ned: number of edges; degs: degree sequence
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_fbnr0(const dg_t *g, double *nr, int method,
    int csepmethod, int *ned, int *degs)
{
  int err, nedges = -1;
  double fb;
  DG_DEFN_(g)

#ifndef DGFB_NEDA
#define DGFB_NEDA 2
#endif
#ifndef DGFB_NEDB
#define DGFB_NEDB (-3)
#endif

  if ( DG_N_ == 1 ) {
    *nr = 1;
    return 1;
  } else if ( DG_N_ == 2 ) {
    if ( g->c[0] ) return -(*nr = 1);
    else return (*nr = 0);
  }

  if ( ned == NULL ) ned = &nedges;

  /* 1. special cases: very loose and very dense graphs */
  if ( !(method & DGSC_NOSPEC) ) {
    fb = dg_fbnr_spec0(g, nr, ned, degs, &err);
    if ( err == 0 ) return fb;
  }

  /* 2. detect clique separators */
  if ( csepmethod != DGCSEP_NULLMETHOD && dg_csep0(g, csepmethod) != 0 ) {
    /* if there is a clique separator, fb = 0
     * directly compute the ring content if needed */
    if (nr != NULL) *nr = dgring_nr0(g, ned, degs);
    return 0;
  }

  /* 3. if the graph is not very dense,
   * use the direct method */
  if ( *ned <= DGFB_NEDA*DG_N_ + DGFB_NEDB
    || DG_N_ > RJWNMAX ) {
    /* ring content is guaranteed to be available here */
    return dgsc_fbnr0(g, nr, method, ned, degs);
  }

  /* 4. finally call the RJW method compute fb,
   * and the direct method for nr */
  if (nr != NULL) *nr = dgring_nr0(g, ned, degs);
  return (double) dgrjw_fb(g);
}



/* free all stock objects */
static void dgrjw_free(void)
{
  if (dgrjw_vsbysize_ != NULL) {
    free(dgrjw_vsbysize_);
    dgrjw_vsbysize_ = NULL;
  }
  dgrjw_vsn_ = dgrjw_vsnmax_ = -1;
  if (dgrjw_fbarr_ != NULL) {
    free(dgrjw_fbarr_);
    dgrjw_fbarr_ = NULL;
  }
  if (dgrjw_aparr_ != NULL) {
    free(dgrjw_aparr_);
    dgrjw_aparr_ = NULL;
  }
  dgrjw_nmax_ = 0;
}



#endif

