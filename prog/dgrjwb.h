#ifndef DGRJW2_H__
#define DGRJW2_H__
/* compute the overall weight of a configuration
 * by a modified RJW method */



#include "dgrjw.h"

#ifdef RJWBLDBL
typedef long double dgrjwb_fb_t;
#else
typedef double dgrjwb_fb_t;
#endif


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
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
      fq[vs * DG_N_ + r1] = (r1 == z1) ? dgrjw_fq(g->c, vs) : 0;

  /* 1. compute the ranked zeta transform of fq */
  for (i = 0; i < DG_N_; i++) {
    bi = MKBIT(i);
    for (vs = 1; vs < vsmax; vs++) {
      if (vs & bi) continue;
      z1 = bitcount(vs) - 1; /* 1 <= z < n */
      for (r1 = 0; r1 <= z1; r1++) /* 1 <= r <= z */
        fq[(vs | bi) * DG_N_ + r1] += fq[vs * DG_N_ + r1];
    }
  }
  for (id = 0; id < vsmax * DG_N_; id++) /* make a copy */
    fc[id] = fq1[id] = fq[id];

  for (fac = -1, k = 2; k <= DG_N_; k++, fac *= -1) {
    /* 2. compute fq^k, collect by the rank r */
    for (vs = 1; vs < vsmax; vs++) {
      id = vs * DG_N_;
      /* the k-fold graph has at most |vs|*k vertices */
      z1 = bitcount(vs) * k;
      r1 = (z1 < DG_N_) ? z1 - 1 : DG_N_ - 1;
      /* since we only use nonempty subsets, |vs| >= 1,
       * and the k-fold subset has at least k members.
       * we loop downward to enable in-place storage */
      for (; r1 >= k - 1; r1--) {
        /* the old (k-1)-fold subset has at least k - 1 members
         * the index of fq satisfies r-t >= k-1  -->  t1 <= r1-k+1 */
        for (x = 0, t1 = r1-k+1; t1 >= 0; t1--)
          x += fq[id + (r1-1-t1)] * fq1[id + t1];
        fq[id + r1] = x;
        fc[id + r1] += x * fac / k;
      }
      fq[id + r1] = 0; /* clear the r = k - 1 case */
    }
  }

  /* 3. compute the ranked Mobius transform */
  for (i = 0; i < DG_N_; i++)
    for (bi = MKBIT(i), vs = 1; vs < vsmax; vs++) {
      if (vs & bi) continue; /* only to subsets without i */
      for (r1 = 0; r1 < DG_N_; r1++) /* each rank is independent */
        fc[(vs | bi) * DG_N_ + r1] -= fc[vs * DG_N_ + r1];
    }

  /* 4. clear zero entries */
  for (vs = 1; vs < vsmax; vs++)
    for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
      /* rounding to integers only works for a hard-sphere liquid */
      fc[vs * DG_N_ + r1] = (r1 == z1) ? round(fc[vs * DG_N_ + r1]) : 0;
      
}



/* remove the contribution of subgraphs with an articulation point at v */
static void dgrjwb_getfb(const dg_t *g, int v,
    dgrjwb_fb_t *fb, dgrjwb_fb_t *fc, dgrjwb_fb_t *fc1)
{
  int i, k, t1, r1, z1;
  dgword_t id, vs, vsmax, bv, bi;
  dgrjwb_fb_t x, fac;
  DG_DEFN_(g)

  vsmax = MKBIT(DG_N_);
  bv = MKBIT(v);

  /* 1. compute the ranked zeta transform of fc */
  for (i = 0; i < DG_N_; i++) {
    if (i == v) continue;
    bi = MKBIT(i);
    for (vs = 1; vs < vsmax; vs++) {
      /* only transform a vertex set with v, but without i */
      if ((vs & bi) || !(vs & bv) || vs == bv) continue;
      z1 = bitcount(vs) - 1; /* 2 <= z < n */
      for (r1 = 0; r1 <= z1; r1++) /* 1 <= r <= z */
        fc[(vs | bi) * DG_N_ + r1] += fc[vs * DG_N_ + r1];
    }
  }
  for (id = 0; id < vsmax * DG_N_; id++) /* make a copy */
    fb[id] = fc1[id] = fc[id];

  for (fac = -1, k = 2; k < DG_N_; k++, fac *= -1) {
    /* 2. compute fc^k, collect by the rank r */
    for (vs = 1; vs < vsmax; vs++) {
      if ( !(vs & bv) || vs == bv ) continue;
      id = vs * DG_N_;
      /* except v, the subset has at least |vs| - 1 >= 1 members
       * the k-fold version has (|vs| - 1) * k members */
      z1 = (bitcount(vs) - 1) * k;
      r1 = (z1 < DG_N_ - 1) ? z1 : DG_N_ - 1;
      /* except v, the k-fold subset has at least k members */
      for (; r1 >= k; r1--) {
        /* except v, the old subset has at least k - 1 members
         * so r - t >= k - 1  -->  t1 <= r1 - k + 1
         * and new addition has at least 1 member except v
         * so t >= 1, --> t1 > 0 */
        for (x = 0, t1 = r1 - k + 1; t1; t1--)
          x += fc[id + r1-t1] * fc1[id + t1];
        fc[id + r1] = x;
        fb[id + r1] += x * fac / k;
      }
      fc[id + r1] = 0; /* clear the r = k - 1 case */
    }
  }

  /* 3. compute the ranked Mobius transform */
  for (i = 0; i < DG_N_; i++) {
    if (i == v) continue;
    bi = MKBIT(i);
    for (vs = 1; vs < vsmax; vs++) {
      /* only transform a vertex set with v, but without i */
      if ((vs & bi) || !(vs & bv) || vs == bv) continue;
      for (r1 = 0; r1 < DG_N_; r1++)
        fb[(vs | bi) * DG_N_ + r1] -= fb[vs * DG_N_ + r1];
    }
  }

  /* 4. clear zero entries */
  for (vs = 1; vs < vsmax; vs++)
    if (vs & bv)
      for (z1 = bitcount(vs) - 1, r1 = 0; r1 < DG_N_; r1++)
        /* rounding to integers only works for a hard-sphere liquid */
        fb[vs * DG_N_ + r1] = (r1 == z1) ? round(fb[vs * DG_N_ + r1]) : 0;
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

  dgrjwb_getfc(g, dgrjwb_arr_, &dgrjwb_arr_[size], &dgrjwb_arr_[size*2]);
  for (v = 0; v < DG_N_; v++) {
    id = (v + 1) % 2;
    dgrjwb_getfb(g, v, &dgrjwb_arr_[size*id], &dgrjwb_arr_[size*(!id)],
      &dgrjwb_arr_[size*2]);
  }
  return dgrjwb_arr_[size*id + size - 1];
}



#endif

