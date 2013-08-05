#ifndef DGRJW_H__
#define DGRJW_H__
/* compute the overall weight of a configuration
 * by the method of Richard J. Wheatley,
 * PRL 110, 200601 (2013) */
#include "dg.h"
#include "dgcsep.h"
#include "dgsc.h"
#include <time.h>
#include <limits.h>


/* for n <= 14, max |fb(n)| = 12! = 479001600 < |INT_MIN| */
#define FBDIRTY INT_MIN
#define FBINVALID(c) ((c) == FBDIRTY)



/* for a configuration, compute sum of all diagrams
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dg_hsfqrjwlow(const code_t *c, code_t vs)
{
  code_t w, br;

  /* if there is a bond, i.e., r(i, j) < 1, then fq = 1 */
  for (w = vs; w; w ^= br) {
    int r = bitfirstlow(w, &br);
    /* if c[r] share vertices with vs, there is bond
     * the regular partition function allows no clash
     * therefore return zero immediately */
    if (c[r] & vs) return 0;
  }
  return 1;
}



/* for a configuration, compute the sum of connected diagrams by
 * Wheatley's recursion formula
 * `c' is the connectivity matrix,  `vs' is the vertex set */
static int dg_hsfcrjwlow(const code_t *c, code_t vs,
    int * RESTRICT fcarr, int * RESTRICT fqarr)
{
  int fc, fc1, fq2;
  code_t ms, ms1, vs1, vs2, s, b, b1;

  b1 = vs & (-vs);
  if ( (vs ^ b1) == 0 ) return 1; /* only one vertex */
  if ( FBINVALID(fc = fqarr[vs]) )
    fqarr[vs] = fc = dg_hsfqrjwlow(c, vs); /* start with fq */
  ms = vs ^ b1; /* the set of vertices (except the lowest vertex) */
  /* loop over subsets of vs, stops when vs == vs1 */
  for (ms1 = 0; ms1 ^ ms;) {
    vs1 = ms1 | b1; /* add vertex 1 to the set */
    vs2 = vs1 ^ vs; /* the complement set */
    if ( FBINVALID(fq2 = fqarr[vs2]) )
      fqarr[vs2] = fq2 = dg_hsfqrjwlow(c, vs2); /* fq of the complement set */
    if (fq2 != 0) {
      if ( FBINVALID(fc1 = fcarr[vs1]) )
        fcarr[vs1] = fc1 = dg_hsfcrjwlow(c, vs1, fcarr, fqarr); /* recursion */
      fc -= fc1 * fq2;
    }
    /* update the subset `ms1' */
    s = ms ^ ms1; /* find the set of empty bits (unused vertices) */
    b = s & (-s); /* find the lowest empty bit (first unused vertex) */
    ms1 = (ms1 | b) & ~(b - 1); /* set this bit, and clear all lower bits */
  }
  return fc;
}



/* compute the sum of connected diagrams by Wheatley's method
 * stand-alone driver */
INLINE int dg_hsfcrjw(const dg_t *g)
{
  static int nmax, *fcarr, *fqarr;
  int i, n = g->n;

  if (fcarr == NULL) {
    nmax = n;
    xnew(fcarr, 1u << nmax);
    xnew(fqarr, 1u << nmax);
  } else if (n > nmax) {
    nmax = n;
    xrenew(fcarr, 1u << nmax);
    xrenew(fqarr, 1u << nmax);
  }
  for (i = 0; i < (1 << n); i++) {
    fcarr[i] = FBDIRTY;
    fqarr[i] = FBDIRTY;
  }
  return dg_hsfcrjwlow(g->c, (1 << n) - 1, fcarr, fqarr);
}



INLINE int dg_hsfarjwlow(const code_t *c, int n, int v, code_t vs,
    int * RESTRICT faarr, int * RESTRICT fbarr);



/* compute the sum of all connected diagrams without the articulation point
   at vertices lower than `v', `vs' is the vertex set
   if v = 0, it returns fc; if v = n, it returns fb */
INLINE int dg_hsfbrjwlow(const code_t *c, int n, int v, code_t vs,
    int * RESTRICT faarr, int * RESTRICT fbarr)
{
  int fb, id, i;
  code_t r, b, bv = (code_t) 1u << v;

  if ((i = bitcount(vs)) <= 1) return 1;
  else if (i == 2) return (c[ bitfirst(vs) ] & vs) ? -1 : 0;

  /* start with the sum of connected diagrams, the first 2^n numbers of
   * fbarr and faarr are used for saving fcarr and fqarr, respectively */
  if ( FBINVALID(fbarr[vs]) )
    fbarr[vs] = dg_hsfcrjwlow(c, vs, fbarr, faarr);
  fb = fbarr[vs];
  /* remove diagrams with the lowest articulation points at i < v */
  for (r = vs & (bv - 1); r; r ^= b) {
    i = bitfirstlow(r, &b);
    id = ((i + 1) << n) + vs;
    if ( FBINVALID(faarr[id]) )
      faarr[id] = dg_hsfarjwlow(c, n, i, vs, faarr, fbarr);
    fbarr[id] = (fb -= faarr[id]);
  }
  return fb;
}



/* compute the sum of all connected diagrams with the articulation point at v
   and no articulation point at any vertex lower than v */
INLINE int dg_hsfarjwlow(const code_t *c, int n, int v, code_t vs,
    int * RESTRICT faarr, int * RESTRICT fbarr)
{
  int fa = 0;
  code_t ms, ms1, vs1, vs2, s, b, bv = 1u << v, b1, id1, id2;

  b1 = vs & (-vs); /* lowest vertex */
  if ( b1 == bv ) { /* if vertex 1 coincide with `v', find the next lowest */
    b1 = vs ^ bv; /* remove the bit */
    b1 = b1 & (-b1);
  }
  ms = vs ^ (b1 ^ bv); /* remove vertices 1 and `v' from the iteration */
  if ( ms == 0 ) return 0; /* no articulated diagram with only two vertices */
  /* loop over subsets of vs, stops when vs == vs1 */
  for (ms1 = 0; ms1 ^ ms;) {
    vs1 = ms1 | (b1 | bv);
    id1 = ((v + 1) << n) + vs1;
    if ( FBINVALID(fbarr[id1]) )
      fbarr[id1] = dg_hsfbrjwlow(c, n, v + 1, vs1, faarr, fbarr);
    if ( fbarr[id1] != 0 ) {
      vs2 = (vs1 ^ vs) | bv; /* complement set of vs */
      id2 = ((v + 1) << n) + vs2;
      if ( FBINVALID(fbarr[id2]) )
        fbarr[id2] = dg_hsfbrjwlow(c, n, v + 1, vs2, faarr, fbarr);
      if ( FBINVALID(faarr[id2]) )
        faarr[id2] = dg_hsfarjwlow(c, n, v, vs2, faarr, fbarr);
      fa += fbarr[id1] * (fbarr[id2] + faarr[id2]);
    }
    /* update the subset `ms1' */
    s = ms ^ ms1; /* find the set of empty bits (unused vertices) */
    b = s & (-s); /* find the lowest empty bit (first unused vertex) */
    ms1 = (ms1 | b) & ~(b - 1); /* set this bit, and clear all lower bits */
  }
  return fa;
}



/* compute the sum of biconnected diagrams by Wheatley's method
 * This is a low level function and the test of clique separator
 * is not done here. */
INLINE int dg_hsfbrjw(const dg_t *g)
{
  static int nmax, *faarr, *fbarr;
  int i, n = g->n;

  if (fbarr == NULL) {
    nmax = n;
    xnew(faarr, (nmax + 1) << nmax);
    xnew(fbarr, (nmax + 1) << nmax);
  } else if (n > nmax) {
    nmax = n;
    xrenew(faarr, (nmax + 1) << nmax);
    xrenew(fbarr, (nmax + 1) << nmax);
  }
  for (i = 0; i < ((nmax + 1) << nmax); i++)
    faarr[i] = fbarr[i] = FBDIRTY;
  return dg_hsfbrjwlow(g->c, n, n, (1u << n) - 1, faarr, fbarr);
}



#define dg_hsfbmixed(g) dg_hsfbmixed0(g, 0)

/* directly compute the sum of biconnected diagrams by various strategies
 * a mixed strategy of dg_hsfbrjw() and dg_rhsc()
 * nocsep: if the graph has been tested with no clique separator */
INLINE int dg_hsfbmixed0(const dg_t *g, int nocsep)
{
  int ned, n = g->n, sgn;

  if ( !nocsep && dg_cliquesep(g) ) return 0;
  ned = dg_nedges(g);
  sgn = ned % 2 ? -1 : 1;
  return (ned >= 2*n - 2) ? dg_hsfbrjw(g) : dg_rhsc_low0(g, 1) * sgn;
}



/* compute the hard shpere weight `fb' by a lookup table */
INLINE int dg_hsfb_lookuplow(int n, unqid_t id)
{
  static int *fb[DGMAP_NMAX + 1]; /* fb of unique diagrams */

  if (fb[n] == NULL) { /* initialize the look-up table */
    dg_t *g;
    dgmap_t *m = dgmap_ + n;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock();

    dgmap_init(m, n);
    if (fb[n] == NULL) xnew(fb[n], m->ng);

    /* loop over unique diagrams */
    g = dg_open(n);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g, &m->first[k]);
      if ( dg_biconnected(g) ) {
        fb[n][k] = dg_cliquesep(g) ? 0 : dg_hsfbrjw(g);
        cnt++;
        nz += (fb[n][k] != 0);
      } else fb[n][k] = 0;
    }
    dg_close(g);
    printf("n %d, computed hard sphere weights of %d/%d biconnected diagrams, %gs\n",
        n, cnt, nz, 1.*(clock() - t0)/CLOCKS_PER_SEC);
  }
  return fb[ n ][ id ];
}



/* compute the total hard shpere weight by a lookup table */
INLINE int dg_hsfb_lookup(const dg_t *g)
{
  return dg_hsfb_lookuplow(g->n, dg_getmapid(g));
}



#define dg_hsfb(g) dg_hsfb0(g, 0)

/* compute the hard-sphere total weight of a configuration
 * nocsep: if the graph has been tested with no clique separator */
INLINE int dg_hsfb0(dg_t *g, int nocsep)
{
  return (g->n <= DGMAP_NMAX) ? dg_hsfb_lookup(g) : dg_hsfbmixed0(g, nocsep);
}



#endif

