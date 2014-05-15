#ifndef DGRJWB_H__
#define DGRJWB_H__
/* compute the overall weight of a configuration
 * by a modified RJW method */



#include "dgrjw.h"

#ifdef  RJWBQUAD
#define RJWBF128 RJWBQUAD
#endif

#ifdef RJWBF128
#include <quadmath.h>
typedef __float128 dgrjwb_fb_t;
#define ROUND(x) roundq(x)

#elif defined(RJWBLDBL)
typedef long double dgrjwb_fb_t;
#define ROUND(x) roundl(x)

#else
typedef double dgrjwb_fb_t;
#define ROUND(x) round(x)
#endif

//#define RJWBSLOW 1

/* compute the ranked Mobius transform */
static void dgrjwb_mobius(dgrjwb_fb_t *arr, int n)
{
  int i, r1;
  dgword_t bi, vs, vsmax = MKBIT(n);

  for (i = 0; i < n; i++)
    for (bi = MKBIT(i), vs = 1; vs < vsmax; vs++) {
      if (vs & bi) continue; /* only to subsets without i */
      for (r1 = 0; r1 < n; r1++) /* each rank is independent */
        arr[(vs | bi) * n + r1] -= arr[vs * n + r1];
    }
}



/* compute the connected sum fc */
static void dgrjwb_getfc(const dg_t *g,
    dgrjwb_fb_t *fc, dgrjwb_fb_t *fq, dgrjwb_fb_t *fq1)
{
  int i, k, z1, t1, r1;
  dgword_t id, vs, vsmax, bi;
  dgrjwb_fb_t x, fac;
  DG_DEFN_(g)

  vsmax = MKBIT(DG_N_);

  /* 0. compute fq */
  for (vs = 0; vs < vsmax; vs++)
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++) {
      fq[vs * DG_N_ + r1] = (r1 == z1) ? dgrjw_fq(g->c, vs) : 0;
    }

  /* 1. compute the ranked zeta transform of fq */
  for (i = 0; i < DG_N_; i++) {
    bi = MKBIT(i);
    for (vs = 1; vs < vsmax; vs++) {
      if (vs & bi) continue;
      z1 = bitcount(vs) - 1; /* 1 <= z <= n - 1 */
      for (r1 = 0; r1 <= z1; r1++) /* 1 <= r <= z */
        fq[(vs | bi) * DG_N_ + r1] += fq[vs * DG_N_ + r1];
    }
  }
#ifdef RJWBSLOW
  for (vs = 0; vs < vsmax; vs++) {
    for (r1 = 0; r1 < DG_N_ - 1; r1++)
      fq1[vs * DG_N_ + r1] = fq[vs * DG_N_ + r1];
    /* the rank-n component of fq1 is reserved for fc */
    if (vs)
      fq1[vs * DG_N_ + DG_N_ - 1] = fq[vs * DG_N_ + bitcount(vs) - 1];
  }
#else
  for (id = 0; id < vsmax * DG_N_; id++) /* make a copy */
    fc[id] = fq1[id] = fq[id];
#endif

  for (fac = -1, k = 2; k <= DG_N_; k++, fac *= -1) {
    /* 2. compute fq^k, collect by the rank r */
    for (vs = 1; vs < vsmax; vs++) {
      id = vs * DG_N_;
      /* the k-fold graph has at most |vs|*k vertices */
      z1 = bitcount(vs);
      r1 = (z1 * k < DG_N_) ? z1 * k - 1 : DG_N_ - 1;
      z1--;
      /* since we only use nonempty subsets, |vs| >= 1,
       * and the k-fold subset has at least k members.
       * we loop downward to enable in-place storage */
      for (; (t1 = r1 - k + 1) >= 0; r1--) {
        /* the old (k-1)-fold subset has at least k - 1 members
         * the index of fq satisfies r-t >= k-1
         *   -->
         * t1 <= r1-k+1 <= n - k, so the t1 = n - 1 is never used */
        /* fc1[t-1] is zero for t > |vs|
         * this condition appears to be unhelpful */
        //if (t1 > z1) t1 = z1;
        for (x = 0; t1 >= 0; t1--)
          x += fq[id + (r1-1-t1)] * fq1[id + t1];
        fq[id + r1] = x;
#ifndef RJWBSLOW
        fc[id + r1] += x * fac / k;
#endif
      }
      fq[id + r1] = 0; /* clear the r = k - 1 case */
    }
#ifdef RJWBSLOW
    for (id = 0; id < vsmax * DG_N_; id++) fc[id] = fq[id];
    dgrjwb_mobius(fc, DG_N_);
    for (vs = 1; vs < vsmax; vs++)
      fq1[vs * DG_N_ + DG_N_ - 1] += fc[vs * DG_N_ + bitcount(vs) - 1] * fac / k;
#endif
  }

#ifdef RJWBSLOW
  for (vs = 1; vs < vsmax; vs++) /* copy fc */
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
      fc[vs * DG_N_ + r1] = (r1 == z1) ? fq1[vs * DG_N_ + DG_N_ - 1] : 0;
#else
  /* 3. compute the ranked Mobius transform */
  dgrjwb_mobius(fc, DG_N_);

  /* 4. clear zero entries */
  for (vs = 1; vs < vsmax; vs++)
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
      /* rounding to integers only works for a hard-sphere liquid */
      fc[vs * DG_N_ + r1] = (r1 == z1) ? ROUND(fc[vs * DG_N_ + r1]) : 0;
#endif
}



/* modified Mobius transform */
static void dgrjwb_mobiusb(dgrjwb_fb_t *arr, int n, int v)
{
  int i, r1;
  dgword_t bvi, bi, bv = MKBIT(v), vs, vsmax = MKBIT(n);

  for (i = 0; i < n; i++) {
    if (i == v) continue;
    bvi = (bi = MKBIT(i)) ^ bv;
    for (vs = 1; vs < vsmax; vs++)
      if ((vs & bvi) == bvi) /* only transform a vertex set with v and i */
        for (r1 = 1; r1 < n; r1++) /* 2 <= r <= n */
          arr[vs * n + r1] -= arr[(vs ^ bi) * n + r1];
  }
}



/* remove the contribution of subgraphs with an articulation point at v */
static void dgrjwb_getfb(const dg_t *g, int v,
    dgrjwb_fb_t *fb, dgrjwb_fb_t *fc, dgrjwb_fb_t *fc1)
{
  int i, k, t1, r1, z1;
  dgword_t id, vs, vsmax, bv, bi, bvi;
  dgrjwb_fb_t x, fac;
  DG_DEFN_(g)

  vsmax = MKBIT(DG_N_);
  bv = MKBIT(v);

  /* 0. clear single-bit diagrams */
  for (i = 0; i < DG_N_; i++)
    fc[MKBIT(i) * DG_N_] = 0;

  /* 1. compute the ranked zeta transform of fc */
  for (i = 0; i < DG_N_; i++) {
    if (i == v) continue;
    bvi = (bi = MKBIT(i)) ^ bv;
    for (vs = 1; vs < vsmax; vs++)
      if ((vs & bvi) == bvi) { /* only transform a vertex set with v and i */
        z1 = bitcount(vs) - 1; /* 2 <= z <= n - 1 */
        for (r1 = 1; r1 <= z1; r1++) /* 2 <= r <= z <= n - 1 */
          fc[vs * DG_N_ + r1] += fc[(vs ^ bi) * DG_N_ + r1];
      }
  }
#ifdef RJWBSLOW
  for (vs = 0; vs < vsmax; vs++) {
    for (r1 = 0; r1 < DG_N_ - 1; r1++)
      fc1[vs * DG_N_ + r1] = fc[vs * DG_N_ + r1];
    if (vs)
      fc1[vs * DG_N_ + DG_N_ - 1] = fc[vs * DG_N_ + bitcount(vs) - 1];
  }
#else
  for (id = 0; id < vsmax * DG_N_; id++) /* make a copy */
    fb[id] = fc1[id] = fc[id];
#endif

  for (fac = -1, k = 2; k < DG_N_; k++, fac *= -1) {
    /* 2. compute fc^k, collect by the rank r */
    for (vs = 1; vs < vsmax; vs++) {
      if ( !(vs & bv) ) continue;
      id = vs * DG_N_;
      /* except v, the subset has at least |vs| - 1 >= 1 members
       * the k-fold version has k ... (|vs| - 1) * k members */
      z1 = bitcount(vs) - 1;
      r1 = (z1 * k < DG_N_ - 1) ? z1 * k : DG_N_ - 1;
      /* except v, the k-fold subset has at least k members */
      for (; (t1 = r1 - k + 1) > 0; r1--) {
        /* except v, the old subset has at least k - 1 members
         * so r - t >= k - 1   -->   t1 <= r1 - k + 1 */
        /* fc1[t-1] is zero for t > |vs|
         * this condition is not helpful */
        //if (t1 > z1) t1 = z1;
        /* and new addition has at least one member except v
         * so t >= 2   -->   t1 > 0 */
        /* r1 - t1 + 1 <= (k - 1) * (|vs| - 1) + 1
         * --> t1 >= r1 - (k - 1) * z1, this condition is unused */
        for (x = 0; t1; t1--)
          x += fc[id + r1-t1] * fc1[id + t1];
        fc[id + r1] = x;
#ifndef RJWBSLOW
        fb[id + r1] += x * fac / k;
#endif
      }
      fc[id + r1] = 0; /* clear the r = k - 1 case */
    }
#ifdef RJWBSLOW
    for (id = 0; id < vsmax * DG_N_; id++) fb[id] = fc[id];
    dgrjwb_mobiusb(fb, DG_N_, v);
    for (vs = 1; vs < vsmax; vs++)
      if (vs & bv)
        fc1[vs * DG_N_ + DG_N_ - 1] += ROUND(fb[vs * DG_N_ + bitcount(vs) - 1] * fac / k);
#endif
  }

#ifdef RJWBSLOW
  for (vs = 1; vs < vsmax; vs++) /* copy fb */
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
      fb[vs * DG_N_ + r1] = (r1 == z1) ? fc1[vs * DG_N_ + DG_N_ - 1] : 0;
#else
  /* 3. compute the ranked Mobius transform */
  dgrjwb_mobiusb(fb, DG_N_, v);

  /* 4. clear zero entries */
  for (vs = 1; vs < vsmax; vs++)
    if (vs & bv)
      for (z1 = bitcount(vs) - 1, r1 = 1; r1 < DG_N_; r1++)
        /* rounding to integers only works for a hard-sphere liquid */
        fb[vs * DG_N_ + r1] = (r1 == z1) ? ROUND(fb[vs * DG_N_ + r1]) : 0;
  /* 3. compute the ranked Mobius transform */
  dgrjwb_mobius(fc, DG_N_);
#endif
}



static int dgrjwb_nmax_;
static dgrjwb_fb_t *dgrjwb_arr_;
/* every thread needs its own memory to do independent calculation
 * so these variables must be thread private */
#pragma omp threadprivate(dgrjwb_nmax_, dgrjwb_arr_)



/* compute the sum of biconnected subgraphs by the modified Wheatley method */
static dgrjwb_fb_t dgrjwb_fb(const dg_t *g)
{
  size_t size = 0;
  int v, id = 0;
  DG_DEFN_(g)

#if !defined(N) || N < 32
  size = ((size_t) DG_N_ << DG_N_); /* n 2^n */
#endif

  if (dgrjwb_arr_ == NULL) {
    dgrjwb_nmax_ = DG_N_;
    xnew(dgrjwb_arr_, size * 3);
    fprintf(stderr, "%4d: dgrjwb allocated %gMB memory for n %d\n",
      inode, sizeof(dgrjwb_fb_t)*size*3./(1024*1024), dgrjwb_nmax_);
  } else if (DG_N_ > dgrjwb_nmax_) {
    dgrjwb_nmax_ = DG_N_;
    xrenew(dgrjwb_arr_, size * 3);
    fprintf(stderr, "%4d: dgrjwb reallocated %gMB memory for n %d\n",
      inode, sizeof(dgrjwb_fb_t)*size*3./(1024*1024), dgrjwb_nmax_);
  }

  dgrjwb_getfc(g, dgrjwb_arr_, dgrjwb_arr_+size, &dgrjwb_arr_[size*2]);
  for (v = 0; v < DG_N_; v++) {
    id = (v + 1) % 2;
    dgrjwb_getfb(g, v, &dgrjwb_arr_[size*id], &dgrjwb_arr_[size*(!id)],
      &dgrjwb_arr_[size*2]);
  }
  return dgrjwb_arr_[size*id + size - 1];
}



#endif

