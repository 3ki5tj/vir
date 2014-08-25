#ifndef DGRING_H__
#define DGRING_H__


/* ring content */



#include "dg.h"
#include "dgutil.h"



#define dgring_perm(g) dgring_perm0(g, 0)


/* directly count the number of ring subgraphs */
static double dgring_perm0(const dg_t *g, int root)
{
  int st[DG_NMAX], top, sttop;
  dgvs_t unused, ccp, ccp0, croot, c, masktop, vstmp;
  dgword_t bi;
  DGVS_DEFIQ_(iq)
  double cnt = 0;
  DG_DEFN_(g)

  if (DG_N_ <= 2) {
    if (DG_N_ <= 1) return 1;
    else return (double) (int) (DGVS_FIRSTWORD(g->c[1]) & 0x1);
  }
  st[0] = root;
  st[1] = -1;
  top = 1;
  /* unused = DG_MASKN_ ^ MKBIT(root); */
  DGVS_MKBITSMASK(unused, DG_N_)
  DGVS_REMOVE(unused, root)
  DGVS_CPY(croot, g->c[root]) /* vertices adjacent to the `root' */

  DGVS_AND2(ccp, unused, croot); /* ccp: set of unused vertices adjacent to st[top-1] */
  sttop = st[top];

  while (1) {
    DGVS_CPY(ccp0, ccp) /* backup ccp in case pushing failed */
    /* construct a set of vertices that satisfy
     * 1. in ccp: unused && connected to st[top-1]
     * 2. indices > sttop
     * the second condition is used to avoid repeated search */
    /* c = ccp & MKINVBITSMASK(sttop + 1); */
    DGVS_MKBITSMASK(masktop, sttop + 1)
    DGVS_MINUS2(c, ccp, masktop) /* exclude vertices in masktop */
    if ( dgvs_nonzero(c) ) {
      DGVS_FIRSTLOW(sttop, c, bi, iq) /* choose `sttop' (or bit `bi') from `c' */
      DGVS_XOR1(unused, bi, iq) /* remove b from the unused vertex set */
      DGVS_AND2(ccp, unused, g->c[sttop]) /* ccp = unused & g->c[sttop]; */
      if (ccp) { /* there are still unused vertices */
        if (top < DG_N_ - 2) { /* push */
          st[top++] = sttop; /* save the current top */
          sttop = -1; /* clear the next level */
          continue;
        }
        /* stay on the highest level, fall through
         * ccp should be a single bit representing the last the vertex */
        /* check ccp and the root are adjacent */
        cnt += dgvs_nonzero( dgvs_and2(vstmp, croot, ccp) );
        // cnt += ((croot & ccp) != 0);
      }
      /* stay on this level */
      DGVS_CPY(ccp, ccp0) /* recover the old ccp */
      DGVS_XOR1(unused, bi, iq) /* add b back to the unused vertices */
      //unused ^= b;
    } else { /* exhausted this level, pop  */
      if (--top == 0) break;
      sttop = st[top];
      DGVS_FLIP(unused, sttop) /* add `b' back to the unused vertices */
      // unused ^= MKBIT(sttop);
      DGVS_AND2(ccp, unused, g->c[st[top - 1]]); /* reconstruct ccp */
    }
  }
  return cnt / 2; /* clockwise & counterclockwise */
}



/* compute the ring content by using combinatorical techniques
 * on the complementary graph */
static double dgring_inv(const dg_t *g)
{
  static int degs[DG_NMAX];
  static int chain[DG_NMAX], nb[DG_NMAX][2];
  static int ed[DG_NMAX * (DG_NMAX - 1) / 2][2];
  static int st[DG_NMAX * (DG_NMAX - 1) / 2 + 2];
#pragma omp threadprivate(degs, chain, nb, ed, st)
  int n = g->n, i, j, k, nj, ied, ned, ked, nch, kch;
  int top, degi, degj, nfree, chi, chj;
  double *fact = dg_facts_, *pow2 = dg_pow2s_;
  double nr = 0, x;

  DG_CALC_FACTS()
  DG_CALC_POW2S()
  /* collect edges of the inverse graph */
  for (ned = 0, i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++)
      if ( !dg_linked(g, i, j) )
        ed[ned][0] = i, ed[ned][1] = j, ned++;
  /* initialize the chain ids and neighbor list */
  for (i = 0; i < n; i++) {
    chain[i] = 0;
    nb[i][0] = nb[i][1] = -1;
  }
  ked = 0;
  nch = 0; /* number of chains */
  kch = 0; /* chain id, ever increasing */
  nr = .5 * fact[n - 1];
  st[ top = 0 ] = -1;
  while ( top >= 0 ) {
    /* search over remaining edges */
    for ( ied = st[top] + 1; ied < ned; ied++ ) {
      if ( (degi = degs[ i = ed[ied][0] ]) >= 2
        || (degj = degs[ j = ed[ied][1] ]) >= 2 )
        continue;
      nfree = (degi == 0) + (degj == 0);
      /* ring detection, if both ends belong to the same chain
       * 1. if the number of chains >= 2, a sub-ring is found
       * 2. if the number of chains == 1, and # of vertices have not
       *    covered all vertices, a sub-ring if also found.
       * We only allow a ring that covers all vertices
       * the number vertices involved is equal to the number of edges
       * plus the number of chains  */
      if ( nfree == 0 && chain[i] == chain[j]
        && (nch >= 2 || nch + ked < n) )
          continue;
      ked++;
      degs[i] = degi + 1;
      degs[j] = degj + 1;
      /* if nb[i][0] != -1 (occupied), we update nb[i][1], otherwise, nb[i][0] */
      nb[i][ (nb[i][0] != -1) ] = j;
      nb[j][ (nb[j][0] != -1) ] = i;
      /* if both vertices are of degree zero, we create a chain
       * if both vertices are of degree one, we join two chains
       * if either is 0 or 1, one chain is extended */
      if (nfree == 2) {
        ++kch;
        chain[i] = kch;
        chain[j] = kch;
        nch++;
      } else if (nfree == 0) { /* joining two existing chains */
        if ( (chi = chain[i]) != (chj = chain[j]) ) {
          //for (k = 0; k < n; k++) if (chain[k] == chj) chain[k] = chi;
          chain[k = j] = chi;
          j = nb[j][ (nb[j][0] == -1) ];
          do {
            chain[j] = chi;
            nj = nb[j][ (nb[j][0] == k) ];
            k = j;
            j = nj;
          } while (j != -1);
        }
        nch--;
      } else if (degi == 0) {
        chain[i] = chain[j];
      } else {
        chain[j] = chain[i];
      }
      st[top] = ied;
      x = (nch == 0) ? 1. : fact[n - ked - 1] * pow2[nch - 1];
      nr += (1 - ked % 2 * 2) * x;
      //for (k = 0; k < n; k++) if ( chain[k] ) printf("%d: %d; ", k, chain[k]); printf("\nx %g, nr %g, edge %d (%d, %d), degi %d, degj %d, chaini %d, chainj %d, nfree %d, nch %d, ked %d, nv %d\n", x, nr, ied, i, j, degi, degj, chain[i], chain[j], nfree, nch, ked, nch + ked); //getchar();
      break;
    }
    if ( ied < ned ) { /* try to push */
      if ( ++top < ned ) {
        /* edges in the stack are arranged in ascending order, so the edge
         * on the next position must have a greater index than this edge */
        st[top] = ied;
        continue;
      }
    }
    --top;
    if ( top >= 0 && (ied = st[top]) >= 0 ) {
      i = ed[ied][0];
      j = ed[ied][1];
      //die_if (nb[i][0] != j && nb[i][1] != j || nb[j][0] != i && nb[j][1] != i, "removing: bad neighbor i %d: %d, %d, j %d: %d, %d\n", i, nb[i][0], nb[i][1], j, nb[j][0], nb[j][1]);
      /* if nb[i][0] != j, we empty nb[i][1], otherwise nb[i][0] */
      nb[i][ (nb[i][0] != j) ] = -1;
      nb[j][ (nb[j][0] != i) ] = -1;
      degi = degs[i] - 1; degs[i] = degi;
      degj = degs[j] - 1; degs[j] = degj;
      //die_if (degi < 0 || degj < 0, "bad degrees %d, %d\n", i, j);
      nfree = (degi == 0) + (degj == 0);
      if ( nfree == 0 ) { /* break the chain */
        /* give vertices connected to j a new chain id */
        k = i;
        ++nch;
        ++kch;
        do {
          chain[j] = kch;
          /* if nb[j][0] == k, choose nb[j][1], otherwise nb[j][0]
           * as the next vertex */
          nj = nb[j][ (nb[j][0] == k) ];
          k = j;
          j = nj;
        } while (nj >= 0);
      } else {
        /* if nfree == 2, removing (i, j) deletes a chain */
        if ( nfree == 2 ) nch--;
        if ( degi == 0 ) chain[i] = 0;
        if ( degj == 0 ) chain[j] = 0;
      }
      ked--;
    }
  }
  return nr;
}



/* the ring content
 * because of the memory limit, we must have n < 32 */
static double dgring_dplow(const dg_t *g, double *arr, unsigned *vsbysize)
{
  int n1, i, j, root;
  dgword_t vs, vs1, vsmax, ms, ms1, id, b, bj;
  double nr;
  DG_DEFN_(g)

  root = DG_N_ - 1;
  n1 = DG_N_ - 1;
  vsmax = (dgword_t) 1 << n1;

  /* zero vertex vs = 0 */
  arr[root] = 1;

  /* one-vertex subset */
  for (id = 1; id <= (dgword_t) DG_N_ - 1; id++) {
    vs = vsbysize[id];
    i = BIT2ID(vs);
    arr[vs * n1 + i] = dg_linked(g, root, i);
  }

  /* general case */
  for ( ; id < vsmax; id++ ) {
    /* loop over vertices in vs */
    for ( vs1 = vs = vsbysize[id]; vs1 != 0; ) {
      vs1 ^= (b = vs1 & (-vs1));
      i = BIT2ID(b);
      ms = vs ^ b;
      /* count paths to i passing all vertices in ms */
      ms1 = DGVS_FIRSTWORD(g->c[i]) & ms;
      for ( nr = 0; ms1 != 0; ) {
        ms1 ^= (bj = ms1 & (-ms1));
        j = BIT2ID(bj);
        nr += arr[ms * n1 + j];
      }
      arr[vs * n1 + i] = nr;
    }
  }

  /* sum over the last link */
  for ( nr = 0, vs1 = DGVS_FIRSTWORD(g->c[root]); vs1 != 0; ) {
    vs1 ^= (b = vs1 & (-vs1));
    i = BIT2ID(b);
    nr += arr[(vsmax - 1) * n1 + i];
  }
  return nr * .5;
}



/* subsets of n vertices sorted by the size */
static unsigned *dgring_vsbysize_;
static int dgring_vsn_ = -1, dgring_vsnmax_ = -1;
#pragma omp threadprivate(dgring_vsbysize_, dgring_vsn_, dgring_vsnmax_)

static double *dgring_nrarr_;
static int dgring_nmax_;
#pragma omp threadprivate(dgring_nrarr_, dgring_nmax_)

/* compute the ring content by dynamic programming */
static double dgring_dp(const dg_t *g)
{
  size_t size = 0;
  DG_DEFN_(g)

#if !defined(N) || N < 32
  size = ((size_t) 1 << (DG_N_ - 1)) * DG_N_;
#endif

  /* allocate memory */
  if (dgring_nrarr_ == NULL) {
    dgring_nmax_ = DG_N_;
    xnew(dgring_nrarr_, size);
  } else if (DG_N_ > dgring_nmax_) {
    dgring_nmax_ = DG_N_;
    xrenew(dgring_nrarr_, size);
  }

  dgring_vsbysize_ = dg_prep_vsbysize(DG_N_ - 1, dgring_vsbysize_,
      &dgring_vsn_, &dgring_vsnmax_);

  return dgring_dplow(g, dgring_nrarr_, dgring_vsbysize_);
}



/* compute the ring content by Karp's algorithm
 * Richard K. Karp, Operations Research Letters
 * Vol. 1, No. 2, April 1982
 * we can safely assume n <= 32, otherwise way too slow */
static double dgring_karp(const dg_t *g)
{
  DG_DEFN_(g)
  int s = DG_N_ - 1, t, r, k;
  dgword_t bs, bt, br;
  dgword_t ms, ms1, ms2, vs;
  double nr, nr1, nr2, ps[2][DG_NMAX] = {{0}};

  bs = (unsigned) 1 << s;
  for ( nr = 0, ms = 1; ms < bs; ms++ ) {
    /* compute the number of walks of length n
     * passing through s and the vertices in ms */
    for (t = 0; t < DG_N_; t++) ps[0][t] = (t == s);
    vs = ms | bs;
    for ( k = 1; k < DG_N_; k++ ) {
      for ( ms1 = vs; ms1 != 0; ) { /* loop over the t */
        BITFIRSTLOW(t, ms1, bt)
        ms1 ^= bt;
        ms2 = DGVS_FIRSTWORD(g->c[t]) & vs;
        for ( nr2 = 0; ms2 != 0; ) { /* loop over r */
          BITFIRSTLOW(r, ms2, br)
          ms2 ^= br;
          nr2 += ps[(k - 1) % 2][r];
        }
        ps[k % 2][t] = nr2;
      }
    }
    ms1 = DGVS_FIRSTWORD(g->c[s]) & ms;
    for ( nr1 = 0; ms1 != 0; ) { /* loop over t */
      BITFIRSTLOW(t, ms1, bt)
      ms1 ^= bt;
      nr1 += ps[(DG_N_ - 1) % 2][t];
    }
    nr += nr1 * (1 - (DG_N_ - bitcount(vs)) % 2 * 2);
  }
  return 0.5 * nr;
}



#define dgring_nr(g) dgring_nr0(g, NULL, NULL)

/* best strategy of computing the ring content */
static double dgring_nr0(const dg_t *g, int *ned, int *degs)
{
  int method = 0;
  DG_DEFN_(g)

  DG_CALC_DEGS(g, ned, degs, dg_nedges_, dg_degs_)

#ifndef DGNR_NEDA
#define DGNR_NEDA 2.5
#endif
#ifndef DGNR_NEDB
#define DGNR_NEDB -0.51
#endif

#ifndef DGNR_NEDC
#define DGNR_NEDC 1.5
#endif
#ifndef DGNR_NEDD
#define DGNR_NEDD -3.01
#endif

#ifndef DGNR_DPNMAX
#define DGNR_DPNMAX 24
#endif

  if ( DG_N_ <= DGNR_DPNMAX ) {
    if (*ned >= DG_N_*(DG_N_ - 1)/2 - (DG_N_*DGNR_NEDC + DGNR_NEDD) ) {
      method = 2;
    } else if (*ned >= DG_N_*DGNR_NEDA + DGNR_NEDB) {
      method = 1;
    } else method = 0;
  } else {
    method = (*ned < DG_N_*(0.25*DG_N_ + 1.) - 5.)  ? 0 : 2;
  }

  if (method == 0) return dgring_perm(g);
  else if (method == 1) return dgring_dp(g);
  else if (method == 2) return dgring_inv(g);
  else return dgring_karp(g); // unused currently
}



/* free all stock objects */
static void dgring_free(void)
{
  if (dgring_vsbysize_ != NULL) {
    free(dgring_vsbysize_);
    dgring_vsbysize_ = NULL;
  }
  dgring_vsn_ = dgring_vsnmax_ = -1;

  if (dgring_nrarr_ != NULL) {
    free(dgring_nrarr_);
    dgring_nrarr_ = NULL;
  }
  dgring_nmax_ = -1;
}



#endif /* DGRING_H__ */

