#ifndef DGCRYR_H__
#define DGCRYR_H__
/* compute the direct correlation function
 * and the cavity function */
#include "dg.h"
#include "dgrjw.h"



/* add or remove the edge between a and b
 * return 0 if the value of fb is zero */
INLINE int dgyr_preproc(dg_t *g, int a, int b, int type, int *eab)
{
  DG_DEFN_(g)

  if ( type == 0 ) { /* compute y(r) */
    /* join a and b, if they are not adjacent */
    if ( !dg_linked(g, a, b) ) {
      dg_link(g, a, b);
      *eab = 1;
    } else {
      *eab = 0;
    }
  } else { /* compute y(r) - t(r) */
    if ( dg_linked(g, a, b) ) {
      dg_unlink(g, a, b);
      if ( !dg_biconnected(g)
        || (type == 2 && !dg_connectedvs(g, MKBITSMASK(DG_N_) ^ MKBIT(a) ^ MKBIT(b))) ) { /* already disqualified */
        dg_link(g, a, b);
        return 0;
      }
      *eab = 1;
    } else {
      *eab = 0;
    }
  }
  return 1;
}



/* remove or add back the edge between a and b,
 * which was added or removed in dgyr_preproc() */
INLINE int dgyr_postproc(dg_t *g, int a, int b, int type, int eab)
{
  if ( type == 0 ) { /* y(r) */
    if ( eab ) dg_unlink(g, a, b); /* disconnect a-b */
    return -1; /* remove the bond between a and b */
  } else { /* y(r) - t(r) or ln y(r) */
    if ( eab ) dg_link(g, a, b);
    return 1;
  }
}


#define dgsc_yriter(g, a, b, type) dgsc_yriter0(g, a, b, type, 0)

/* iteratively compute the signed star content for
 * y(r) [type = 0] or y(r) - t(r) [type = 1] or ln y(r) [type = 2]
 *
 * 0. y(r) is computed as follows
 * First, join the root vertices a-b, if they are not adjacent in g
 * Then, the signed star content of a-b is equal to
 * the number of biconnected subgraphs of g
 * with even number of missing edges (except a-b)
 * less the number of biconnected subgraphs of g
 * with odd number of missing edges (except a-b)
 *
 * 1. y(r) - t(r) = c(r) - f(r) y(r) is computed as follows
 * First, disjoin the root vertices a and b, if they are adjacent in g
 * The remaining is the same.
 *
 * 2. Similar to case 1, but also exclude graphs with {a, b} being
 * articulation pair.
 *
 * On return g should be the same
 * This function is adapted from dgsc_iter0.h from dgsc.h
 * */
INLINE double dgsc_yriter0(dg_t *g, int a, int b, int type, int noproc)
{
  static int ed[DG_NMAX * (DG_NMAX - 1)/2][2]; /* edges */
  static int st[DG_NMAX * (DG_NMAX - 1)/2 + 2]; /* state */
#pragma omp threadprivate(ed, st)
  int n = g->n, top, ied, med, vi, vj, degi, degj, eab;
  dgword_t vsnab = MKBITSMASK(n) ^ MKBIT(a) ^ MKBIT(b);
  double fb;
  int ned = -1, degs[DG_NMAX];

  /* make sure a < b */
  if ( a > b ) vi = a, a = b, b = vi;

  die_if (b < 0, "(%d-%d) is not a valid edge\n", a, b);

  if (!noproc && dgyr_preproc(g, a, b, type, &eab) == 0) return 0;

  /* collect edges to be switched on or off */
  ned = dg_degs(g, degs);
  for (med = 0, vi = 0; vi < n - 1; vi++)
    if ( degs[vi] > 2 )
      for (vj = vi + 1; vj < n; vj++)
        if ( degs[vj] > 2 && dg_linked(g, vi, vj)
             && !(vi == a && vj == b) ) { /* the edge a-b is fixed */
          ed[med][0] = vi;
          ed[med][1] = vj;
          med++;
        }

  st[ top = 0 ] = -1;
  fb = 1 - ned % 2 * 2;
  while (top >= 0) {
    /* search over the remaining edges */
    for (ied = st[top] + 1; ied < med; ied++) {
      /* removing an edge adjacent to a degree-2 vertex
       * makes a graph not biconnected */
      if ( (degi = degs[ vi = ed[ied][0] ]) < 3
        || (degj = degs[ vj = ed[ied][1] ]) < 3 )
        continue;
      dg_unlink(g, vi, vj);
      if ( dg_biconnected(g)
          && (type != 2 || dg_connectedvs(g, vsnab)) ) {
        /* if we want the bridge diagrams, va and vb cannot be a separation pair */
        degs[vi] = degi - 1, degs[vj] = degj - 1;
        st[top] = ied;
        --ned;
        fb += 1 - ned % 2 * 2;
        break;
      }
      dg_link(g, vi, vj);
    }
    if (ied < med) { /* try to push */
      if (++top < med) {
        /* edges in the stack are arranged in ascending order, so the edge
         * on the next position must have a greater index than this edge */
        st[top] = ied; /* clean up the next level */
        continue;
      } /* otherwise, fall through and pop */
    }
    --top; /* pop */
    if (top >= 0 && (ied = st[top]) >= 0) { /* erase the edge */
      vi = ed[ied][0]; vj = ed[ied][1];
      dg_link(g, vi, vj);
      degs[vi]++, degs[vj]++;
      ned++;
    }
  }

  return noproc ? fb : fb * dgyr_postproc(g, a, b, type, eab);
}




/* compute, for y(r), the modified Boltzmann weight,
   the sum of all Mayer diagrams with an f-bond between a and b
  `c' is the connectivity matrix, `vs' is the vertex set */
INLINE int dgrjw_yrfq(dgvs_t *c, dgword_t vs, int a, int b)
{
  dgword_t w, bi;

  /* if there is a bond, i.e., r(i, j) < 1, then fq = 0 */
  for (w = vs; w; w ^= bi) {
    int i = bitfirstlow(w, &bi);
    /* if c[i] share vertices with vs, there is bond
     * the Boltzmann weight = \prod_(ij) e_ij, and it allows no clash
     * therefore return zero immediately
     * exclude, however, the edge a-b */
    if ( i == a && (DGVS_FIRSTWORD(c[i]) & (vs & ~b)) ) return 0;
    else if ( i == b && (DGVS_FIRSTWORD(c[i]) & (vs & ~a))) return 0;
    else if (DGVS_FIRSTWORD(c[i]) & vs) return 0;
  }
  return 1; /* (vs != 0); */
}



/* compute, for y(r), the modified Boltzmann weight (fq) and
 * the sum of connected diagrams (fc)
 * also a few values for fa and fb
 * `fqarr' extends to `faarr', `fcarr' extends to `fbarr'
 * `vsbysize[i]' gives the ith subset of `g->n' vertices,
 *    where the subsets are sorted by the number of vertices contained
 * `aparr[vs]' is the index array of the lowest articulation point + 1
 *    in the diagram induced by the vertex set `vs'
 * by Wheatley's recursion formula, for all diagrams */
static void dgrjw_yrprepare(const dg_t *g, int a, int b, int type,
    dgrjw_fb_t *fcarr, dgrjw_fb_t *fqarr,
    unsigned *vsbysize, dgrjw_ap_t *aparr)
{
  dgrjw_fb_t fc;
  dgword_t vsab, vs, vs1, ms, ms1, ms2, b1, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */

  vsab = MKBIT(a) | MKBIT(b);

  /* we loop from smaller subsets to larger ones */
  /* diagrams with zero or one vertex */
  for (id = 0; id <= (unsigned) DG_N_; id++) {
    vs = vsbysize[id];
    die_if (bitcount(vs) > 1, "n %d, id %d, bad count(%#x) = %d\n", DG_N_, id, vs, bitcount(vs));
    fcarr[vs] = 1;
    fqarr[vs] = 1;
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

    /* start with the modified Boltzmann weight */
    if ( (vs & vsab) == vsab )
      fc = dgrjw_yrfq(g->c, vs, a, b);
    else /* if the vertex set involves neither a nor b */
      fc = dgrjw_fq(g->c, vs);
    fqarr[vs] = fc;

    if ( !dg_connectedvs(g, DGVS_W2VS(vs)) ) {
      /* disconnected subset */
      fcarr[vs] = 0;
      aparr[vs] = 0; /* fb(vs) == 0 */
    } else { /* connected diagram */
      if ( (vs & vsab) == vsab && type == 0 ) { /* vs contains a and b */
        ms = vs ^ vsab; /* exclude the roots from the vertex set */
        /* loop over proper subsets of `vs'
         * the first component contains `b1'
         * ms1 is connected, contains a and b
         * ms2 is arbitrary */
        for (ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
          fc -= fcarr[ms1 ^ vsab] * fqarr[ms2];
          DGVS_INCa(ms1, ms2); /* update the subset `ms1' */
        }
      } else { /* vs does not contain a and b */
        b1 = vs & (-vs); /* lowest vertex */
        ms = vs ^ b1;
        for (ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
          fc -= fcarr[ms1 ^ b1] * fqarr[ms2];
          DGVS_INCa(ms1, ms2); /* update the subset `ms1' */
        }
      }
      fcarr[vs] = fc;

      /* check articulation points */
      aparr[vs] = (dgrjw_ap_t) (DG_N_ + 1);
      for ( vs1 = vs; vs1 != 0; ) {
        b1 = vs1 & (-vs1);
        vs1 ^= b1;
        if ( !fcarr[vs ^ b1] ) {
          /* fb(vs, v) == 0 for v = BIT2ID(b1) + 1 to n - 1 */
          aparr[vs] = (dgrjw_ap_t) (BIT2ID(b1) + 1);
          break;
        }
      }
    } /* connected / disconnected diagrams */
  }
}



/* return fb, assuming dgrjw_yrprepare() has been called */
static dgrjw_fb_t dgrjw_yriter(const dg_t *g, int a, int b, int type,
    dgrjw_fb_t *fbarr, int offset, unsigned *vsbysize, dgrjw_ap_t *aparr)
{
  int v, iold, inew;
  size_t idold, idnew, jdold, jdnew;
  dgrjw_fb_t fa;
  dgword_t vs, vsab, vsnew, ms, ms1, ms2, b1, bv, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */
  vsab = MKBIT(a) | MKBIT(b);

  /* one- and two-vertex subsets, biconnected == connected */
  for ( id = 0; id <= (dgword_t) (DG_N_ * (DG_N_ + 1) / 2); id++ ) {
    vs = vsbysize[id];
    fbarr[ vsmax | vs ] = fbarr[ vs ];
  }

  /* loop over the position of the first clique separator */
  for ( v = 0; v < DG_N_; v++ ) {
    iold = (v + offset) % 2;
    inew = (v + 1 + offset) % 2;
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
          if ( (vs & vsab) == vsab && type == 0 ) {
            /* in this case, vs1 must contain both a and b */
            ms = vs & ~vsab & ~bv; /* the looper set excludes a, b, v */
            jdnew = idnew | vsab | bv;
          } else {
            /* in this case, vs1 must contain the lowest vertex v1 except v */
            ms = vs & ~bv;
            b1 = ms & (-ms); /* lowest vertex */
            ms &= ~b1; /* remove the fixed vertices `b1' and `bv' from `vs' */
            jdnew = idnew | b1 | bv;
          }
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
  return fbarr[vsmax * ((DG_N_ + offset) % 2) + vsmax - 1];
}



/* remove graphs with {a, b} being a separation pair */
static dgrjw_fb_t dgrjw_yrrmseppair(const dg_t *g, int a, int b,
    dgrjw_fb_t *fbarr, dgrjw_fb_t *fdarr, unsigned *vsbysize)
{
  dgrjw_fb_t fa;
  dgword_t vs, vsab, ms, ms1, ms2, b1, id, vsmax;
  DG_DEFN_(g)

  vsmax = 1u << DG_N_; /* 2^n */
  vsab = MKBIT(a) | MKBIT(b);

  /* three-vertex subsets */
  for ( id = DG_N_ * (DG_N_ + 1) / 2 + 1 ; id < vsmax; id++ ) {
    vs = vsbysize[id];
    fa = 0;
    /* skip if the vertex set does not contain v */
    if ( (vsab & vs) == vsab ) {
      /* in this case, vs1 must contain the lowest vertex v1 except v */
      ms = vs ^ vsab;
      b1 = ms & (-ms); /* lowest vertex */
      ms &= ~b1; /* remove the fixed vertices `b1' and `bv' from `vs' */
      for ( ms1 = 0; (ms2 = ms1 ^ ms) != 0; ) {
        fa += fdarr[(vsab ^ b1) ^ ms1] * fbarr[vsab ^ ms2];
        /* update the subset `ms1' */
        DGVS_INCa(ms1, ms2);
      }
    }
    fdarr[vs] = fbarr[vs] - fa; /* set fb */
  }
  return fdarr[vsmax - 1];
}



#define dgrjw_yrfb(g, a, b, type) dgrjw_yrfb0(g, a, b, type, 0)

/* compute the sum of biconnected diagrams by Wheatley's method
 * This is a low level function and the test of clique separator
 * is not done here.
 * The memory requirement is 2^(n-20)*21MB */
static dgrjw_fb_t dgrjw_yrfb0(dg_t *g, int a, int b, int type, int noproc)
{
  DG_DEFN_(g)
  size_t size = 0;
  int eab, offset = 0;
  dgrjw_fb_t fb;

#if !defined(N) || N < 32  /* this condition should always be true */
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

  /* add/remove the edge between a and b */
  if ( !noproc && dgyr_preproc(g, a, b, type, &eab) == 0 ) return 0;
  /* pre-compute all fc and fq values */
  dgrjw_yrprepare(g, a, b, type, dgrjw_fbarr_, dgrjw_fbarr_ + size,
                   dgrjw_vsbysize_, dgrjw_aparr_);
  /* remove graphs in which a and b form an articulation pair */
  if ( type == 2 ) {
    dgrjw_yrrmseppair(g, a, b, dgrjw_fbarr_, dgrjw_fbarr_ + size,
                      dgrjw_vsbysize_);
    offset = 1;
  }
  /* remove articulated graphs */
  fb = dgrjw_yriter(g, a, b, type, dgrjw_fbarr_, offset,
                    dgrjw_vsbysize_, dgrjw_aparr_);
  /* remove/add back the edge between a and b */
  return noproc ? fb : fb * dgyr_postproc(g, a, b, type, eab);
}



#if DGVS_ONEWORD



/* find clique separator of two vertices except a b */
static dgvsref_t dg_csep2xp(const dg_t *g, int a, int b)
{
  int i;
  dgvs_t ci, maski, mask;
  dgword_t bi, bj;
  DG_DEFN_(g)
  DG_DEFMASKN_()

  if ( a > b ) i = a, a = b, b = i;
  /* loop for the first vertex */
  for (i = 1; i < DG_N_; i++) {
    mask = MKBITSMASK(i);
    if ( i == b ) mask ^= MKBIT(a);
    bi = MKBIT(i);
    maski = DG_MASKN_ ^ bi;

    /* loop over vertices connected to `i', and with indices less than `i' */
    for (ci = g->c[i] & mask; ci; ci ^= bj) {
      bj = ci & (-ci);
      /* if the cluster is not connected, we have found a clique separator */
      if ( !dg_connectedvs(g, maski ^ bj) )
        return bi ^ bj; /* bi | bj */
    }
  }
  return 0;
}



#else /* !DGVS_ONEWORD */



/* find clique separator of two vertices except {a, b} */
static dgvsref_t dg_csep2xp(const dg_t *g, int a, int b)
{
  int i;
  dgvs_t ci, maski, maskci, mask;
  dgword_t bi, bj;
  DGVS_DEFIQ_(iq)
  DGVS_DEFIQ_(jq)
  DG_DEFN_(g)
  DG_DEFMASKN_()

  if ( a > b ) i = a, a = b, b = i;
  DGVS_CPY(maski, DG_MASKN_)
  DGVS_CLEAR(maskci)
  maskci[0] = 1; /* all vertices lower than i */

  /* loop for the first vertex */
  for (i = 1; i < DG_N_; i++) {
    DGVS_MKBIT(i, bi, iq) /* bi = MKBIT(i); */
    DGVS_XOR1(maski, bi, iq)
    /* maskci: vertices connected to `i', and with indices less than `i' */
    DGVS_AND2(ci, g->c[i], maskci)
    if (i == b) { DGVS_REMOVE(ci, a) }
    DGVS_OR1(maskci, bi, iq) /* update maskci */

    while ( dgvs_nonzero(ci) ) {
      DGVS_FIRSTBIT(ci, bj, jq)
      DGVS_XOR1(ci, bj, jq) /* remove `bj' from `ci' */
      DGVS_XOR1(maski, bj, jq) /* `maski' now lacks the `j' bit */

      /* if the cluster is not connected, we have found a clique separator */
      if ( !dg_connectedvs(g, maski) ) {
        static dgvs_t csep;
        #pragma omp threadprivate(csep)
        DGVS_CLEAR(csep)
        DGVS_OR1(csep, bi, iq)
        DGVS_OR1(csep, bj, jq)
        return csep;
      }
      DGVS_XOR1(maski, bj, jq) /* add back the `j' bit */
    }
    DGVS_XOR1(maski, bi, iq) /* add back the `i' bit */
  }
  return 0;
}


#endif /* DGVS_ONEWORD */




#define dg_yrfb(g, a, b, type) dg_yrfb0(g, a, b, type, NULL, NULL)

/* compute y(r), or y(r) - t(r), or log y(r) */
INLINE double dg_yrfb0(dg_t *g, int a, int b, int type,
                       int *ned, int *degs)
{
  double fb = 0;
  int eab, nedges = -1;
  DG_DEFN_(g)

  if ( ned == NULL ) ned = &nedges;

  if (*ned < 0) {
    DG_CALC_DEGS(g, ned, degs, dg_nedges_, dg_degs_); /* prepare degrees */
  }

  if ( dgyr_preproc(g, a, b, type, &eab) == 0 ) return 0;

  if ( dg_csep2xp(g, a, b) != 0 ) {
    fb = 0;
    goto END;
  }

  /* if the graph is not very dense, use the direct method */
  if ( *ned <= 2*DG_N_ - 2 || DG_N_ > RJWNMAX ) {
    fb = dgsc_yriter0(g, a, b, type, 1);
  } else {
    fb =  (double) dgrjw_yrfb0(g, a, b, type, 1);
  }

END:
  return fb * dgyr_postproc(g, a, b, type, eab);
}



/* compute y(r), or y(r) - t(r), or log y(r) */
INLINE double *dg_yrfball(double *fb, dg_t *g, int type)
{
  int a, b, id, ned, degs[DG_NMAX];
  DG_DEFN_(g)

  ned = dg_degs(g, degs); /* prepare degrees */
  id = 0;
  for ( a = 0; a < DG_N_; a++ )
    for ( b = a + 1; b < DG_N_; b++, id++ )
      fb[id] = dg_yrfb0(g, a, b, type, &ned, degs);
  return fb;
}


/* histogram for correlation functions */
typedef struct {
  double xmin;
  double xmax;
  double dx;
  double sphr;
  double vol;
  double *arr;
  int dim;
  int order;
  int npt; /* number of points along r */
  unsigned long nsampi;
  double nsampd;
} hscr_t;



INLINE hscr_t *hscr_open(real xmax, real dx, int dim, int order)
{
  hscr_t *hs;
  int i;

  xnew(hs, 1);
  hs->xmin = 0;
  hs->xmax = xmax;
  hs->dx = dx;
  hs->dim = dim;
  hs->order = order;
  hs->npt = (int)(xmax/dx + 0.99999999);
  hs->nsampi = 0;
  hs->nsampd = 0;
  hs->sphr = (dim % 2 + 1) * dim;
  for (i = 2 + dim % 2; i <= dim; i += 2 )
    hs->sphr *= M_PI * 2 / i;
  hs->vol = hs->sphr / dim;
  xnew(hs->arr, hs->npt);
  return hs;
}



/* increase the count of the histogram */
INLINE void hscr_inc(hscr_t *hs, int cnt)
{
  if ( (hs->nsampi += cnt) > 1000000000 ) {
    hs->nsampd += (double) hs->nsampi;
    hs->nsampi = 0;
  }
}



/* add an entry, but do not update the sample size */
INLINE void hscr_add0(hscr_t *hs, double x, double fb)
{
  int ix = (int)((x - hs->xmin)/hs->dx);
  if (ix < hs->npt) hs->arr[ix] += fb;
}



/* add an entry and update the sample size */
INLINE void hscr_add(hscr_t *hs, double x, double fb)
{
  hscr_add0(hs, x, fb);
  hscr_inc(hs, 1);
}



INLINE int hscr_save(hscr_t *hs, const char *fn, double norm)
{
  FILE *fp;
  int i, imax;
  double y, r;

  /* Zn/V^(n-1) --> Zn */
  norm *= pow(hs->vol, hs->order - 1);
  /* divide it by (n - 2)! for permutations of the black dots */
  for ( i = 2; i <= hs->order - 2; i++ ) norm /= i;
  for ( imax = hs->npt - 1; imax >= 0; imax-- )
    if ( fabs(hs->arr[imax]) > 0 ) break;
  if ( imax < 0 ) return -1;
  /* flush the counter */
  hs->nsampd += (double) hs->nsampi;
  hs->nsampi = 0;
  xfopen(fp, fn, "w", return -1);
  fprintf(fp, "# %d %d %d %.0f %g %g %g %g\n",
      hs->dim, hs->order, imax + 1,
      hs->nsampd, hs->dx, norm, hs->sphr, hs->vol);
  for ( i = 0; i <= imax; i++ ) {
    r = (i + .5) * hs->dx;
    y = hs->arr[i]/(hs->nsampd * hs->sphr * pow(r, hs->dim - 1) * hs->dx);
    if (norm > 0) y *= norm;
    fprintf(fp, "%12.7f %22.14e\n", r, y);
  }
  fclose(fp);
  return 0;
}



INLINE void hscr_close(hscr_t *hs)
{
  if (!hs) return;
  if (hs->arr) free(hs->arr);
  memset(hs, 0, sizeof(*hs));
}



#endif /* defined DGCRYR_H__ */

