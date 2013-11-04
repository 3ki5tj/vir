#ifndef DGAUT_H__
#define DGAUT_H__
/* canonical label */

#include "dg.h"


/* set parameteters for the program `nauty' (No AUTomorphism, Yes?) */
#define WORDSIZE        DG_WORDBITS
#define MAXN            DG_NMAX
#define ONE_WORD_SETS 1 /* use one-word set when possible */
#include "nau0s.h"



/* canonical label */
INLINE dg_t *dg_canlabel(dg_t *gout, const dg_t *gin)
{
  static dgword_t g0[DG_NMAX], gc[DG_NMAX];
  static int lab[DG_NMAX], ptn[DG_NMAX], orbits[DG_NMAX];
  static DEFAULTOPTIONS_GRAPH(options);
  static statsblk stats;
#pragma omp threadprivate(g0, gc, lab, ptn, orbits, options, stats)
  int i;
  DG_DEFN_(gin);

  /* nauty uses the highest bit for the first index */
  for (i = 0; i < DG_N_; i++) g0[i] = bitreverse(gin->c[i]);
  options.getcanon = TRUE;
  densenauty(g0, lab, ptn, orbits, &options, &stats, 1, DG_N_, gc);

  gout->n = DG_N_;
  for (i = 0; i < DG_N_; i++) gout->c[i] = bitreverse(gc[i]);
  return gout;
}



typedef struct {
  int n; /* number of vertices */
  int nc; /* number of sets in the partition */
  dgword_t cs[DG_NMAX]; /* cs[i] is the ith vertex set */
  int cnt[DG_NMAX]; /* cnt[i] is the number of vertices in set i */
} dgpart_t; /* an ordered partition, a sequence of vertex sets */



INLINE void dgpart_unit(dgpart_t *part, int n)
{
  part->n = n;
  part->nc = 1;
  part->cs[0] = MKBITSMASK(n);
  part->cnt[0] = n;
}



#define dgpart_print(p) dgpart_fprint(p, stdout)
#define dgpart_errprint(p) dgpart_fprint(p, stderr)

/* print out the partition */
INLINE void dgpart_fprint(const dgpart_t *part, FILE *fp)
{
  int ic, k;
  dgword_t vs, b;

  for (ic = 0; ic < part->nc; ic++) {
    for (vs = part->cs[ic]; vs; vs ^= b) {
      BITFIRSTLOW(k, vs, b);
      fprintf(fp, "%d ", k);
    }
    fprintf(fp, "| ");
  }
  fprintf(fp, "\n");
}



/* check if a partition part is valid
 * set `equitable' to 1 check if the partition is equitable */
INLINE int dgpart_check(const dgpart_t *part, int equitable, const dg_t *g)
{
  int ic, sz;
  dgword_t all = 0, vs;

  for (ic = 0; ic < part->nc; ic++) {
    vs = part->cs[ic];
    if (vs & all) {
      fprintf(stderr, "cell %d has a previous element %d\n",
          ic, BIT2ID(vs & all));
      dgpart_errprint(part);
      return 1;
    }
    all |= vs;
    sz = bitcount(vs);
    if (sz != part->cnt[ic] || sz == 0) {
      fprintf(stderr, "bad cell %d size %d, or it mismatches %d\n",
          ic, sz, part->cnt[ic]);
      dgpart_errprint(part);
      return 2;
    }
  }
  if (all != MKBITSMASK(part->n)) {
    fprintf(stderr, "partition does not cover all vertices\n");
    dgpart_errprint(part);
    return 3;
  }

  if (equitable) { /* check if the diagram is equitable */
    /* loop over the shatterers */
    for (ic = 0; ic < part->nc; ic++) {
      int jc, deg0, deg, k, k0;
      dgword_t vsi, vsj, b;
      vsi = part->cs[ic];
      /* loop over the shatterees */
      for (jc = 0; jc < part->nc; jc++) {
        if (part->cnt[jc] == 1) continue;
        deg0 = k0 = -1;
        /* loop over vertices in the cell jc */
        for (vsj = part->cs[jc]; vsj; vsj ^= b) {
          BITFIRSTLOW(k, vsj, b);
          deg = dg_degvs(g, k, vsi);
          if (deg0 < 0) {
            deg0 = deg;
            k0 = k;
          } else if (deg != deg0) {
            fprintf(stderr, "cell %d (point %d) is not equitable w.r.t. cell %d, degs %d vs %d (point %d)\n",
                jc, k, ic, deg, deg0, k0);
            dgpart_print(part);
            dg_print(g);
            return 99;
          }
        }
      }
    }
  }
  return 0;
}



/* return an equitable partition of the graph that is compatible with part
 * An equitable partition is a sequence of vertex subsets or cells {Vi},
 * such that for any Vi and Vj (i can be equal to j) and
 * for any vk and vl in Vi, deg(vk, Vj) == deg(vk, Vi) */
INLINE void dg_equipart(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g);
  int i, i0, imax, j, k, deg, ni;
  dgword_t vsj, vsi, vs, b;
  dgword_t vswdeg[DG_NMAX]; /* vswdeg[i] is the vertex set with degree i */

  for (j = 0; j < part->nc; j++) { /* loop over shatterers */
    vsj = part->cs[j];
    //printf("shatter j %d\n", j); dgpart_print(part); getchar();

    /* loop over shatterees
     * In this process, the shatterer j is treated as fixed
     * so once the cell i is shattered properly, its subcells
     * cannot be further shattered by j, this means that we
     * only need to loop to the current number of cells imax
     * even if part->nc grows latter on.
     * Note, the shatterer j may be changed during the process.
     * However, if Vj_old cannot shatter a cell i*, but the
     * new Vj_new (a subset of Vj_old) can shatter this cell i*
     * then the subset Vj_old - Vj_new can shatter this cell i* as well
     * And Vj_old - Vj_new contains all newly created subcells appended
     * at the end of the list.
     * Thus, it is okay to jump to the next j even if j is changed */
    imax = part->nc; /* the current number of cells */
    i0 = -1; /* the first shattered cell */
    for (i = 0; i < imax; i++) {
      vsi = part->cs[i];
      ni = part->cnt[i]; /* number of vertices */
      //die_if (ni == 0, "set %d is empty\n", i);
      if (ni == 1) { /* one-vertex cell, nothing to shatter */
        continue;
      } else {
        dgword_t nzdegs, nz, bdeg;
        /* we use count sort to shatter vsi, that is,
         * we sort cells according to their degrees with respect to vsj
         * then each nonempty bags labeled by the degree gives the shatterer */
        /* loop over vertices in the shatteree */
        nzdegs = 0;
        for (vs = vsi; vs; vs ^= b) {
          BITFIRSTLOW(k, vs, b);
          deg = dg_degvs(g, k, vsj);
          bdeg = MKBIT(deg);
          if (bdeg & nzdegs) { /* vswdeg[deg] already has something */
            vswdeg[deg] |= b;
          } else {
            vswdeg[deg] = b;
            nzdegs |= bdeg; /* mark degree `deg' has been occupied */
          }
        }
        /* if there is only one occupied cell, no shattering */
        if (bitcount(nzdegs) == 1) continue;
        /* now nonzero vswdeg[x] collects a subset of vsi with deg x
         * we now loop over occupied degrees */
        k = 0; /* number of shattered cells */
        for (nz = nzdegs; nz; nz ^= bdeg) {
          BITFIRSTLOW(deg, nz, bdeg);
          if (k == 0) { /* replace cell i by the first partition */
            part->cnt[i] = bitcount( part->cs[i] = vswdeg[deg] );
            k = part->nc; /* the next set begins here */
          } else {
            part->cnt[k] = bitcount( part->cs[k] = vswdeg[deg] );
            part->nc = ++k;
            if (i0 < 0) i0 = i; /* update i0, if cell i has been shattered */
          }
        }
      }
      if (part->nc == DG_N_) return;
    } /* end of the loop over shatterees */
  } /* end of the loop over shatterers */
  //die_if (dgpart_check(part, 1, g) != 0, "equipart made a bad partition %d", DG_N_);
}



#define dg_printperm(perm, n) dg_fprintperm(stdout, perm, n)

/* print a permutation */
INLINE void dg_fprintperm(FILE *fp, const int *perm, int n)
{
  int i;

  for (i = 0; i < n; i++)
    fprintf(fp, "%d ", perm[i]);
  fprintf(fp, "\n");
}



/* check if a permuation of n is good, return 0 on success */
INLINE int dg_checkperm(const int *perm, int n)
{
  dgword_t c = 0, b;
  int i;

  for (i = 0; i < n; i++) {
    if (perm[i] < 0 || perm[i] >= n) {
      dg_fprintperm(stderr, perm, n);
      return 1;
    }
    b = MKBIT(perm[i]);
    if (c & b) {
      fprintf(stderr, "permutation has two %d\n", perm[i]);
      dg_fprintperm(stderr, perm, n);
      return 2;
    }
    c |= b;
  }
  return 0;
}



/* get a graph that is compatible with the equitable partition
 * if recur is true, try to artificially break the equitable partion
 *   until it becomes discrete
 * the speed of using recur or not is very similar, recommended */
INLINE void dg_permequipart(int *perm, const dg_t *g, int recur)
{
  DG_DEFN_(g);
  dgpart_t part;
  int ic, k, i;
  dgword_t vs, b;

  /* initially set part as the unit partition */
  dgpart_unit(&part, DG_N_);

  if (recur) {
    /* recursively call dg_equipart() until part becomes discrete */
    while (part.nc < DG_N_) {
      dg_equipart(&part, g);
      if (part.nc == DG_N_) break;
      /* find the first unfixed cell */
      for (ic = 0; ic < part.nc; ic++)
        if (part.cnt[ic] > 1)
          break;
      /* artificially break the first cell */
      vs = part.cs[ic];
      b = vs & (-vs);
      part.cs[part.nc] = vs ^ b;
      part.cnt[part.nc] = part.cnt[ic] - 1;
      part.cs[ic] = b;
      part.cnt[ic] = 1;
      part.nc += 1;
    }
    for (ic = 0; ic < DG_N_; ic++) {
      perm[ BIT2ID(part.cs[ic]) ] = ic;
    }
  } else { /* single pass */
    dg_equipart(&part, g);
    /* loop over cells in the partition */
    for (i = 0, ic = 0; ic < part.nc; ic++) {
      for (vs = part.cs[ic]; vs; vs ^= b) {
        BITFIRSTLOW(k, vs, b);
        perm[ k ] = i++;
      } /* loop over vertices in a partition */
    } /* loop over partitions */
  }
}



/* get a permutation compatible with the degree sequence */
INLINE void dg_permdegseq(int *perm, const dg_t *g)
{
  DG_DEFN_(g);
  int i, k, deg;
  dgword_t cvs[DG_NMAX], vs, b;
#define DG_DEGMIN 0 /* can be 2 for biconnected graphs */

  for (i = DG_DEGMIN; i < DG_N_; i++) cvs[i] = 0;
  /* use count sort to collect vertices with the same degrees */
  for (b = 1, i = 0; i < DG_N_; i++, b <<= 1)
    cvs[dg_deg(g, i)] |= b;
  i = 0;
  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    for (vs = cvs[deg]; vs; vs ^= b) { /* loop over this partition */
      BITFIRSTLOW(k, vs, b);
      perm[ k ] = i++;
    }
  }
}



/* get a partition in which each cell contains all vertices of the same degree */
INLINE void dg_partdegseq(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g);
  int i, deg, ic;
  dgword_t cvs[DG_NMAX], vs, b;
#define DG_DEGMIN 0 /* can be 2 for biconnected graphs */

  for (i = DG_DEGMIN; i < DG_N_; i++) cvs[i] = 0;
  /* use count sort to collect vertices with the same degrees */
  for (b = 1, i = 0; i < DG_N_; i++, b <<= 1)
    cvs[dg_deg(g, i)] |= b;
  /* construct a partition accordingly */
  part->n = g->n;
  ic = 0;
  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    if ((vs = cvs[deg]) == 0) continue;
    part->cs[ic] = vs;
    part->cnt[ic] = bitcount(vs);
    ic++;
  }
  part->nc = ic;
}



/* rearrange vertices in `gin' according to the permuation `perm' */
INLINE dg_t *dg_permutate(dg_t * RESTRICT gout, const dg_t * RESTRICT gin, const int *perm)
{
  DG_DEFN_(gin);
  int i, j;
  dgword_t ci, ci_out;

  gout->n = DG_N_;
  dg_empty(gout); /* adding this appears to make the function faster? */
  for (i = 0; i < DG_N_; i++) {
    for (ci = gin->c[ i ], ci_out = 0, j = 0; j < DG_N_; j++, ci >>= 1)
      if (ci & 0x1)
        ci_out |= MKBIT(perm[j]);
    gout->c[ perm[i] ] = ci_out;
  }
  return gout;
}



/* check if a is a vertex permutation of b */
INLINE int dg_checkiso(const dg_t *a, const dg_t *b, const int *perm)
{
  DG_DEFN_(b);
  int i, j;
  dgword_t cai, cbi;

  for (i = 0; i < DG_N_; i++)
    for (cai = a->c[i], cbi = b->c[perm[i]], j = 0; j < DG_N_; j++)
      if (((cai >> j) & 0x1) != ((cbi >> perm[j]) & 0x1))
        return 1;
  return 0;
}



/* get a representative isomorphic graph `gout' of the input graph `gin'
 * `level' represent the depth of the search
 * if level == 3, it returns the unique canonical label of the graph
 *   all isomorphic graphs share the same canonical label
 * if level == 2, it returns the first node from McKay's search of the
 *    canonical label
 * if level == 1, it returns a graph that is compatible with equitable
 *    partition (hierarchical degree sequence) of the graph
 * if level == 0, it returns a graph compatible with the degree sequence
 *
 * In all cases, `gout' is isomorphic to `gin' since the former is produced
 *  from some vertex permuation of the latter.
 * But isomorphic graphs may not be mapped to the same `gout' if level < 3
 * with a higher `level' the redundance reduces, but the speed goes down.
 * gout->n will be proper set here
 * */
INLINE dg_t *dg_repiso(dg_t *gout, const dg_t *gin, int level)
{
  int perm[DG_NMAX + 1];

  if (level == 0) {
    DG_DEFN_(gin);
    gout->n = DG_N_;
    return dg_copy(gout, gin);
  }

  if (level >= 4 || level < 0) { /* return the canonical label */
    return dg_canlabel(gout, gin);
  }

  if (level == 1) { /* degree sequence */
    dg_permdegseq(perm, gin);
  } else if (level == 2) { /* equitable partition */
    dg_permequipart(perm, gin, 0);
  } else if (level == 3) { /* equitable partition */
    dg_permequipart(perm, gin, 1);
  }
  //die_if (dg_checkperm(perm, gin->n) != 0, "bad permutation level %d\n", level);
  dg_permutate(gout, gin, perm);
  //die_if (dg_checkiso(gout, gin, perm) != 0, "bad isomorphism!\n");
  return gout;
}



/* enumerate all permutations compatible with a partition  */
INLINE int dg_part2permls(int *perms, int nmax, const dgpart_t *part)
{
  int i, ic, k, pcnt, n = part->n, nc = part->nc;
  dgword_t vs, vs1, b, unused;
  int st[DG_NMAX], stdefid[DG_NMAX], top = 0;
  dgword_t stvs[DG_NMAX], stavail[DG_NMAX];

  /* vertices available in each stack position */
  for (i = 0, ic = 0; ic < nc; ic++) {
    int defid = (part->cnt[ic] == 1) ? bitfirst(part->cs[ic]) : -1;
    for (k = 0; k < part->cnt[ ic ]; k++) {
      stvs[i] = part->cs[ ic ];
      stdefid[i] = defid;
      i++;
    }
  }
  //dgpart_print(part);
  //for (i = 0; i < n; i++) printf("%d %d\n", i, stdefid[i]);
  //die_if (i != n, "bad partition i %d != n %d\n", i, n);

  pcnt = 0;
  top = 0;
  unused = mkbitsmask(n);
  stavail[top] = stvs[top];
  while (1) {
    /* determine the vertex on the top */
    vs = stavail[top];
    if (vs) { /* we still have vertices on this level */
      BITFIRSTLOW(k, vs, b);
      st[top] = k;
      //printf("push top %d, k %d, stavail %s, unused %s\n", top, k, dg_strset1(stavail[top]), dg_strset1(unused) );
      stavail[top] = vs ^ b; /* remove b from the unused list */
      if (top < n - 2) { /* push */
        unused ^= b;
        if ((vs1 = stvs[top + 1]) == stvs[top]) { /* same cell */
          stavail[top + 1] = unused & vs1; /* copy the unused list */
        } else { /* different cell */
          stavail[top + 1] = vs1;
        }
        top++;
        continue;
      }
      /* here top == n - 2, this means that we are at the next to
       * the top level, then no need to go up further */
      /* need to determine the last vertex */
      if (unused ^ b) { /* the last vertex is in the same cell */
        st[top + 1] = BIT2ID(unused ^ b);
      } else { /* otherwise the last cell must be a discrete cell */
        st[top + 1] = stdefid[top + 1];
      }
      //dgpart_print(part); dg_printperm(st, n); getchar();
      /* TODO: remove this */
      //die_if (dg_checkperm(st, n) != 0, "bad permutation n %d\n", n);
      /* save the permutation */
      for (k = 0; k < n; k++)
        perms[pcnt * n + st[k]] = k;
      if (++pcnt >= nmax) return pcnt;
      /* do not increase top, so that we can explore other choices
       * on this level */
    } else { /* no vertex left on this level, pop */
      if (--top < 0) break;
      unused |= MKBIT(st[top]);
      //printf("pop top %d, k %d, stavail %s, unused %s\n", top, st[top], dg_strset1(stavail[top]), dg_strset1(unused) );
    }
  }
  return pcnt;
}



static int *dgaut_perms_, dgaut_nperms_; /* permutation */
#pragma omp threadprivate(dgaut_perms_, dgaut_nperms_)

/* generate a list of connectivity codes of graphs isomorphic to g
 * corresponding to the isomorphism level `level' */
INLINE int dg_repisocodels(dgword_t *codes, int cwords, int nmax,
    const dg_t *g, int level, dg_t *ng)
{
  dgpart_t part;
  int ncnt, i, ic, j, n = g->n;

  if (level == 0) { /* unit partition, exchange all vertices */
    dgpart_unit(&part, n);
  } else if (level == 1) {
    dg_partdegseq(&part, g);
  } else if (level == 2) {
    dg_equipart(&part, g);
  } else {
    die_if (level < 0 || level > 2, "do not support (%d)\n", level);
  }

  /* allocate spaces for permutations */
  if (n * nmax > dgaut_nperms_) {
    dgaut_nperms_ = n * nmax;
    if (dgaut_perms_) free(dgaut_perms_);
    xnew(dgaut_perms_, dgaut_nperms_);
  }
  /* obtain a list of possible permutations */
  ncnt = dg_part2permls(dgaut_perms_, nmax, &part);
  /* convert each permutation to a code */
  for (ic = 0, i = 0; i < ncnt; i++) {
    //printf("i %d, ", i); dg_printperm(dgaut_perms_ + i * n, n);
    dg_permutate(ng, g, dgaut_perms_ + i * n);
    //dg_print(gout);
    dg_encode(ng, codes + ic * cwords);
    /* avoid duplicated codes see if the code is new */
    for (j = 0; j < ic; j++)
      if (dgcode_eq(codes + ic * cwords, codes + j * cwords, cwords))
        break;
    /* increment ic only if the code is new */
    if (j >= ic) ic++;
  }
  //getchar();
  return ic;
}



#endif /* DGAUT_H__ */

