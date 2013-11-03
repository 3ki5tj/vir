#ifndef DGRJW_H__
#define DGRJW_H__
/* compute the overall weight of a configuration
 * by the method of Richard J. Wheatley, PRL 110, 200601 (2013) */
#include "dgcsep.h"
#include "dgsc.h"
#include <time.h>
#include <limits.h>


#if defined(RJW32) || (defined(N) && N <= 14 && DG_WORDBITS == 32)

typedef int32_t dgrjw_fb_t;
/* for n <= 14, max |fb(n)| = 12! = 479001600 < |FBDIRTY| */
#ifndef RJWNMAX
#define RJWNMAX 14
#endif
#define FBDIRTY ((dgrjw_fb_t) 0x80808080) /* -2139062144 */
#define FBINVALID(c) ((c) == FBDIRTY)

#else /* 64-bit RJW */

typedef int64_t dgrjw_fb_t;
/* for n <= 22, max |fb(n)| = 20! = 2.4e18 < |FBDIRTY| */
/* Note: the threshold is based on the limit 20! = 2.4e18
 *  and 21! cannot be contained in a 64-bint integer
 * The code may however fail before that due to the cancellation
 *  in intermediate steps
 * Another limit is that memory > 2^(n + 1) * (n + 1) */
#ifndef RJWNMAX
#define RJWNMAX 22
#endif
#define FBDIRTY ((dgrjw_fb_t) 0x8080808080808080ull) /* -9187201950435737472 */
#define FBINVALID(c) ((c) == FBDIRTY)
#endif


/* for a configuration, compute sum of all diagrams
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dg_hsfq_rjwlow(const dgword_t *c, dgword_t vs)
{
  dgword_t w, b;

  /* if there is a bond, i.e., r(i, j) < 1, then fq = 0 */
  for (w = vs; w; w ^= b) {
    /* if c[i] share vertices with vs, there is bond
     * the Boltzmann weight = \prod_(ij) e_ij, and it allows no clash
     * therefore return zero immediately */
    if (c[bitfirstlow(w, &b)] & vs) return 0;
  }
  return 1;
}



#define DGRJW_INC(ms1, ms2) { dgword_t nms2_, lbit_, hbits_; \
  nms2_ = (dgword_t) (-(ms2)); /* -ms2 and ms2 share bits higher than the lowest bit of ms2 */ \
  lbit_ = (ms2) & nms2_; /* the & yields the lowest bit lowbit_ */ \
  hbits_ = (ms2) ^ nms2_; /* the collections of bits higher than lbit_ */ \
  (ms1) &= hbits_; /* this wipes out bits lower or equal to lbit_ */ \
  (ms1) |= lbit_;  /* this add the lowest bit */ }

/* here is an alternative versions, which is slower */
#define DGRJW_INCa(ms1, ms2) { dgword_t lbit_; \
    lbit_ = (ms2) & (-ms2); /* find the lowest empty bit (first unused variable vertex) */ \
    (ms1) ^= lbit_; /* add this bit */ \
    (ms1) &= ~(b - 1); /* clear all lower bits */ }

/* below is another alternative version, which is even slower */
#define DGRJW_INCb(ms1, ms2) { int i = BITFIRSTNZ(ms2); \
  ms1 >>= i; /* clear the lower bits by shifting i bits */ \
  ms1 ^= 1; /* add the lowest bit, `^' can be replaced by `+' or `|' with little difference */ \
  ms1 <<= i; /* shift back */ }



/* for a configuration, compute the sum of connected diagrams by
 * Wheatley's recursion formula
 * `c' is the connectivity matrix,  `vs' is the vertex set */
INLINE dgrjw_fb_t dg_hsfc_rjwlow(const dgword_t *c, dgword_t vs,
    dgrjw_fb_t * RESTRICT fcarr, dgrjw_fb_t * RESTRICT fqarr)
{
  dgrjw_fb_t fc, fc1, fq2;
  dgword_t ms1, ms2, b1;

  b1 = vs & (-vs);
  if ( (vs ^ b1) == 0 ) return 1; /* only one vertex */
  if ( FBINVALID(fc = fqarr[vs]) )
    fqarr[vs] = fc = dg_hsfq_rjwlow(c, vs); /* start with fq */
  vs ^= b1; /* remove vertex b1 from vs, vs = vs - {b1} */
  /* before the above statement, `vs' is the set of all vertices
   * since `b1' is fixed, after excluding `b1' in the statement,
   * `vs' is the set of *variable* vertices
   * `ms1' is the subset of variable vertices
   * `ms2' is the complement set of variable vertices */
  /* loop over subsets ms1 of vs, stops when the complement set ms2 == 0 */
  for (ms1 = 0; (ms2 = ms1 ^ vs) != 0; ) {
    if ( FBINVALID(fq2 = fqarr[ms2]) )
      fqarr[ms2] = fq2 = dg_hsfq_rjwlow(c, ms2); /* fq of the complement set */
    if (fq2 != 0) {
      dgword_t vs1 = ms1 | b1; /* vs1 = ms1 + {b1} */
      if ( FBINVALID(fc1 = fcarr[vs1]) )
        fcarr[vs1] = fc1 = dg_hsfc_rjwlow(c, vs1, fcarr, fqarr); /* recursion */
      fc -= fc1 * fq2;
    }
    /* update the subset `ms1' */
    DGRJW_INC(ms1, ms2);
  }
  return fc;
}



#if 0
/* compute the sum of connected diagrams by Wheatley's method
 * stand-alone driver */
INLINE dgrjw_fb_t dg_hsfc_rjw(const dg_t *g)
{
  static int nmax;
  static dgrjw_fb_t *fcarr, *fqarr;
  int i;
  DG_DEFN_(g);
  DG_DEFMASKN_();

  if (fcarr == NULL) {
    nmax = DG_N_;
    xnew(fcarr, 1u << (nmax + 1));
    fqarr = fcarr + (1u << nmax);
  } else if (DG_N_ > nmax) {
    nmax = DG_N_;
    xrenew(fcarr, 1u << (nmax + 1));
    fqarr = fcarr + (1u << nmax);
  }
  for (i = 0; (unsigned) i < (1u << DG_N_); i++) fcarr[i] = FBDIRTY;
  for (i = 0; (unsigned) i < (1u << DG_N_); i++) fqarr[i] = FBDIRTY;
  return dg_hsfc_rjwlow(g->c, DG_MASKN_, fcarr, fqarr);
}
#endif



INLINE dgrjw_fb_t dg_hsfa_rjwlow(const dgword_t *c, int n, int v, dgword_t vs,
    dgrjw_fb_t * RESTRICT faarr, dgrjw_fb_t * RESTRICT fbarr);



/* compute the sum of all connected diagrams without the articulation point
   at vertices lower than `v', `vs' is the vertex set
   if v = 0, it returns fc; if v = n, it returns fb */
INLINE dgrjw_fb_t dg_hsfb_rjwlow(const dgword_t *c, int n, int v, dgword_t vs,
    dgrjw_fb_t * RESTRICT faarr, dgrjw_fb_t * RESTRICT fbarr)
{
  int i;
  dgrjw_fb_t fb, fa;
  dgword_t r, b, bv = MKBIT(v), id;

  if ((i = bitcount(vs)) <= 1) {
    return 1;
  } else if (i == 2) {
    return (c[ BITFIRSTNZ(vs) ] & vs) ? -1 : 0;
  }

  /* start with the sum of connected diagrams, the first 2^n numbers of
   * fbarr and faarr are used for saving fcarr and fqarr, respectively */
  if ( FBINVALID(fb = fbarr[vs]) ) {
    fb = dg_hsfc_rjwlow(c, vs, fbarr, faarr);
    fbarr[vs] = fb;
  }
  /* remove diagrams with the lowest articulation points at i < v */
  for (r = vs & (bv - 1); r; r ^= b) {
    BITFIRSTLOW(i, r, b);
    id = ((i + 1) << DG_N_) | vs; /* (i + 1) * 2^n + vs
                                   * `|' is equivalent to `+' */
    if ( FBINVALID(fa = faarr[id]) ) {
      fa = dg_hsfa_rjwlow(c, DG_N_, i, vs, faarr, fbarr);
      faarr[id] = fa;
    }
    fbarr[id] = (fb -= fa);
  }
  return fb;
}



/* compute the sum of all connected diagrams with the articulation point at v
   and no articulation point at any vertex lower than v */
INLINE dgrjw_fb_t dg_hsfa_rjwlow(const dgword_t *c, int n, int v, dgword_t vs,
    dgrjw_fb_t * RESTRICT faarr, dgrjw_fb_t * RESTRICT fbarr)
{
  dgrjw_fb_t fa = 0, fb, fa2, fb2;
  dgword_t ms1, ms2, vs1, bv = MKBIT(v), b1, b1v, id0, id;

  b1 = vs & (-vs); /* lowest vertex */
  if ( b1 == bv ) { /* if vertex 1 coincide with `v', find the next lowest */
    b1 = vs ^ bv; /* remove `bv' from the vertex set */
    b1 = b1 & (-b1);
  }
  b1v = b1 ^ bv;
  vs ^= b1v; /* remove the fixed vertices `b1' and `bv' from `vs' */
  /* `vs' is the set of *variable* vertices from now on */
  if ( vs == 0 ) return 0; /* no articulated diagram with only two vertices */
  id0 = ((dgword_t) (v + 1) << DG_N_); /* (v + 1) * 2^n */
  /* `id0' is the offset for vertex v */
  /* loop over subsets of vs, stops when vs == vs1 */
  for (ms1 = 0; (ms2 = (ms1 ^ vs)) != 0; ) {
    vs1 = ms1 | b1v; /* add the two fixed vertices */
    id = id0 | vs1; /* `|' is equivalent to + */
    if ( FBINVALID(fb = fbarr[id]) ) { /* compute fb if necessary */
      fb = dg_hsfb_rjwlow(c, DG_N_, v + 1, vs1, faarr, fbarr);
      fbarr[id] = fb;
    }
    if ( fb != 0 ) {
      dgword_t vs2 = ms2 | bv; /* unused variable vertices + the articulation point `bv'
                              * `|' is equivalent to `+' here */
      id = id0 | vs2; /* `|' is equivalent to + */
      fb2 = fbarr[id];
      if ( FBINVALID(fb2 = fbarr[id]) ) {
        fb2 = dg_hsfb_rjwlow(c, DG_N_, v + 1, vs2, faarr, fbarr);
        fbarr[id] = fb2; /* save the fb value */
      }
      if ( FBINVALID(fa2 = faarr[id]) ) {
        fa2 = dg_hsfa_rjwlow(c, DG_N_, v, vs2, faarr, fbarr);
        faarr[id] = fa2; /* save the fa value */
      }
      fa += fb * (fb2 + fa2);
    }
    /* update the subset `ms1' */
    DGRJW_INC(ms1, ms2);
  }
  return fa;
}



/* compute the sum of biconnected diagrams by Wheatley's method
 * This is a low level function and the test of clique separator
 * is not done here. */
INLINE dgrjw_fb_t dg_hsfb_rjw(const dg_t *g)
{
  static int nmax;
  static dgrjw_fb_t *faarr, *fbarr;
/* every thread needs its own memory to do independent calculation
 * so these variables must be thread private */
#pragma omp threadprivate(nmax, faarr, fbarr)

  DG_DEFN_(g);
  DG_DEFMASKN_();
  size_t size;

  /* the memory requirement is 2^(n + 1) * (n + 1) * sizeof(dgrjw_fb_t) */
  if (fbarr == NULL) {
    /* even if we know N, we still use dynamic allocation
     * because this function might not be called at all
     * in high dimensions */
    nmax = DG_N_;
    size = (nmax + 1) << nmax; /* (nmax + 1) * 2^nmax */
    xnew(faarr, size * 2);
    fbarr = faarr + size;
    fprintf(stderr, "%4d: dgrjw allocated %gMB memory for n %d\n", inode,
        size * 2. * sizeof(dgrjw_fb_t) / (1024*1024), nmax);
  }
#ifndef N /* if N is fixed no need to reallocate */
  else if (DG_N_ > nmax) {
    nmax = DG_N_;
    size = (nmax + 1) << nmax; /* (nmax + 1) * 2^nmax */
    xrenew(faarr, size * 2);
    fbarr = faarr + size;
    fprintf(stderr, "%4d: dgrjw reallocated %gMB memory for n %d\n", inode,
        size * 2. * sizeof(dgrjw_fb_t) / (1024*1024), nmax);
  }
#endif
  size = ((size_t) (DG_N_ + 1) << DG_N_); /* (n + 1) * 2^n */
  /* every byte of FBDIRTY is 0x80, so we can use memset() */
  memset(faarr, 0x80, size * sizeof(faarr[0]));
  memset(fbarr, 0x80, size * sizeof(fbarr[0]));
  return dg_hsfb_rjwlow(g->c, DG_N_, DG_N_, DG_MASKN_, faarr, fbarr);
}



#define dg_hsfb_mixed(g) dg_hsfb_mixed0(g, 0, NULL, NULL)

/* directly compute the sum of biconnected diagrams by various strategies
 * a mixed strategy of dg_hsfb_rjw() and dg_rhsc()
 * nocsep = 1 means that the graph has been tested with no clique separator
 *        = 0 means it MAY have clique separators
 * *ned: number of edges; degs: degree sequence
 * if ned != NULL and *ned <= 0, both *ned and degs[] are computed on return */
INLINE double dg_hsfb_mixed0(const dg_t *g,
    int nocsep, int *ned, int *degs)
{
  int err, nedges = -1;
  DG_DEFN_(g);
  double sc;

  if ( ned == NULL ) ned = &nedges;
  sc = dg_rhsc_spec0(g, nocsep, 1, ned, degs, &err);
  if ( err == 0 ) {
    return DG_SC2FB(sc, *ned);
  } else if ( *ned <= 2*DG_N_ - 3 || DG_N_ > RJWNMAX) {
    return DG_SC2FB(dg_rhsc_directlow(g), *ned);
  } else { /* hsfb_rjw() requires 2^(n + 1) * (n + 1) memory */
    return (double) dg_hsfb_rjw(g);
  }
}



#ifdef DGMAP_EXISTS
#define dg_hsfb_lookup(g) dg_hsfb_lookuplow(g->n, dg_getmapid(g))


/* compute the hard shpere weight `fb' by a lookup table
 * the return can always fit into an intger */
INLINE double dg_hsfb_lookuplow(int n, unqid_t id)
{
  static double *fb[DGMAP_NMAX + 1]; /* fb of unique diagrams */
#pragma omp threadprivate(fb)

  if (DG_N_ <= 1) return 1;

  /* initialize the look-up table */
  if (fb[DG_N_] == NULL) {
    dg_t *g;
    dgmap_t *m = dgmap_ + DG_N_;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock();

    dgmap_init(m, DG_N_);
    xnew(fb[DG_N_], m->ng);

    /* loop over unique diagrams */
    g = dg_open(DG_N_);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g, &m->first[k]);
      if ( dg_biconnected(g) ) {
        fb[DG_N_][k] = (double) (dg_cliquesep(g) ? 0 : dg_hsfb_rjw(g));
        cnt++;
        nz += (fabs(fb[DG_N_][k]) > .5);
      } else fb[DG_N_][k] = 0;
    }
    dg_close(g);
    fprintf(stderr, "%4d: n %d, computed hard sphere weights of %d/%d biconnected diagrams, %gs\n",
        inode, DG_N_, cnt, nz, 1.*(clock() - t0)/CLOCKS_PER_SEC);
  }
  return fb[ DG_N_ ][ id ];
}
#endif /* defined(DGMAP_EXISTS) */


#define dg_hsfb(g) dg_hsfb0(g, 0, NULL, NULL)

/* compute the hard-sphere total weight of a configuration
 * nocsep: if the graph has been tested with no clique separator */
INLINE double dg_hsfb0(dg_t *g, int nocsep, int *ned, int *degs)
{
#ifdef DGMAP_EXISTS
  if (g->n <= DGMAP_NMAX)
    return dg_hsfb_lookup(g);
#endif /* defined(DGMAP_EXISTS) */

  return dg_hsfb_mixed0(g, nocsep, ned, degs);
}



#endif

