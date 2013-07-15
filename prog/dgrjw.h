#ifndef DGRJW_H__
#define DGRJW_H__
/* compute the overall weight of a configuration
 * by the method of Richard J. Wheatley,
 * PRL 110, 200601 (2013) */
#include "dg.h"
#include "dgcsep.h"
#include <time.h>
#include <limits.h>



/* for a configuration, compute sum of all diagrams
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dg_hsfqrjwlow(const dg_t *g, code_t vs)
{
  code_t w, br;

  /* if there is a bond, or r(i, j) < 1, then the fq = 1 */
  for (w = vs; w; w ^= br) {
    int r = bitfirstlow(w, &br);
    /* if g->c[r] share vertices with vs, there is bond */
    if (g->c[r] & vs) return 0;
  }
  return 1;
}



/* for a configuration, compute sum of connected diagrams by Wheatley's method
  `c' is the connectivity matrix,  `vs' is the vertex set */
static int dg_hsfcrjwlow(const dg_t *g, code_t vs, int *fcarr)
{
  int fc, fq2;
  code_t ms, ms1, vs1, vs2, s, b, b1;

  b1 = vs & (-vs);
  if ( (vs ^ b1) == 0 ) return 1; /* only one vertex */
  fc = dg_hsfqrjwlow(g, vs); /* start with fq */
  ms = vs ^ b1; /* the set of vertices (except the lowest vertex) */
  /* loop over subsets of vs, stops when vs == vs1 */
  for (ms1 = 0; ms1 ^ ms;) {
    vs1 = ms1 | b1;
    vs2 = vs1 ^ vs;
    fq2 = dg_hsfqrjwlow(g, vs2); /* fq of the complement set */
    if (fq2 != 0) {
      if (fcarr[vs1] == INT_MIN)
        fcarr[vs1] = dg_hsfcrjwlow(g, vs1, fcarr); /* recursion */
      fc -= fcarr[vs1] * fq2;
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
  static int nmax, *fcarr;
  int i, n = g->n;

  if (fcarr == NULL) {
    nmax = n;
    xnew(fcarr, 1u << (nmax - 1));
  } else if (n > nmax) {
    nmax = n;
    xrenew(fcarr, 1u << (nmax - 1));
  }
  for (i = 0; i < (1 << (n - 1)); i++) fcarr[i] = INT_MIN;
  return dg_hsfcrjwlow(g, (1 << n) - 1, fcarr);
}



INLINE int dg_hsfarjwlow(const dg_t *g, int v, code_t vs, int *faarr, int *fbarr);



/* compute the sum of all connected diagrams without the articulation point
   at vertices lower than `v', `vs' is the vertex set
   if v = 0, it returns fc; if v = n, it returns fb */
INLINE int dg_hsfbrjwlow(const dg_t *g, int v, code_t vs, int *faarr, int *fbarr)
{
  int fb, id, i, n = g->n;
  code_t r, b, bv = (code_t) 1u << v;

  if ((i = bitcount(vs)) <= 1) return 1;
  else if (i == 2) return (g->c[ bitfirst(vs) ] & vs) ? -1 : 0;

  /* start with the sum of connected diagrams */
  if (fbarr[vs] == INT_MIN)
    fbarr[vs] = dg_hsfcrjwlow(g, vs, fbarr);
  fb = fbarr[vs];
  /* remove diagrams with the lowest articulation points at i < v */
  for (r = vs & (bv - 1); r; r ^= b) {
    i = bitfirstlow(r, &b);
    id = ((i + 1) << n) + vs;
    if (faarr[id] == INT_MIN)
      faarr[id] = dg_hsfarjwlow(g, i, vs, faarr, fbarr);
    fbarr[id] = (fb -= faarr[id]);
  }
  return fb;
}



/* compute the sum of all connected diagrams with the articulation point at v
   and no articulation point at any vertex lower than v */
INLINE int dg_hsfarjwlow(const dg_t *g, int v, code_t vs, int *faarr, int *fbarr)
{
  int fa = 0, n = g->n;
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
    if ( fbarr[id1] == INT_MIN )
      fbarr[id1] = dg_hsfbrjwlow(g, v + 1, vs1, faarr, fbarr);
    if ( fbarr[id1] != 0 ) {
      vs2 = (vs1 ^ vs) | bv; /* complement set of vs */
      id2 = ((v + 1) << n) + vs2;
      if ( fbarr[id2] == INT_MIN )
        fbarr[id2] = dg_hsfbrjwlow(g, v + 1, vs2, faarr, fbarr);
      if ( faarr[id2] == INT_MIN )
        faarr[id2] = dg_hsfarjwlow(g, v, vs2, faarr, fbarr);
      fa += fbarr[id1] * (fbarr[id2] + faarr[id2]);
    }
    /* update the subset `ms1' */
    s = ms ^ ms1; /* find the set of empty bits (unused vertices) */
    b = s & (-s); /* find the lowest empty bit (first unused vertex) */
    ms1 = (ms1 | b) & ~(b - 1); /* set this bit, and clear all lower bits */
  }
  return fa;
}




/* compute the sum of biconnected diagrams by Wheatley's method */
INLINE int dg_hsfbrjw(const dg_t *g)
{
  static int nmax, *faarr, *fbarr;
  int i, n = g->n;

  if ( dg_cliquesep(g) ) return 0;
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
    faarr[i] = fbarr[i] = INT_MIN;
  return dg_hsfbrjwlow(g, n, (1u << n) - 1, faarr, fbarr);
}



/* compute the total hard shpere weight by a lookup table */
INLINE int dg_hsfb_lookup(const dg_t *g)
{
  static int *fb[DGMAP_NMAX + 1]; /* biconnectivity of unique diagrams */
  int n = g->n;
  dgmap_t *m = dgmap_ + n;
  code_t c;

  if (fb[n] == NULL) { /* initialize the look-up table */
    dg_t *g1;
    int k, cnt = 0, nz = 0;
    clock_t t0 = clock(), t1;

    if (n >= 8) printf("n %d: initializing...\n", n);
    dgmap_init(m, n); /* compute the permutation mapping */
    if (fb[n] == NULL) xnew(fb[n], m->ng);

    t1 = clock();
    if (n >= 8) printf("n %d: diagram-map initialized %gs\n",
        n, 1.*(t1 - t0)/CLOCKS_PER_SEC);
    /* loop over unique diagrams */
    g1 = dg_open(n);
    for (cnt = 0, k = 0; k < m->ng; k++) {
      dg_decode(g1, &m->first[k]);
      if ( dg_biconnected_lookup(g1) ) { /* use the look-up version */
        fb[n][k] = dg_hsfbrjw(g1);
        cnt++;
        nz += (fb[n][k] != 0);
      } else fb[n][k] = 0;
    }
    dg_close(g1);
    printf("n %d, computed hard sphere weights of %d/%d biconnected diagrams, %gs\n",
        n, cnt, nz, 1.*(clock() - t1)/CLOCKS_PER_SEC);
  }
  dg_encode(g, &c);
  return fb[n][ m->map[c] ]; /* m->map[c] is the id of the unique diagram */
}



/* compute the hard-sphere total weight of a configuration */
INLINE int dg_hsfb(dg_t *g)
{
  return (g->n <= DGMAP_NMAX) ? dg_hsfb_lookup(g) : dg_hsfbrjw(g);
}




#endif

