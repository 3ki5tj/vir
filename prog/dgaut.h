#ifndef DGAUT_H__
#define DGAUT_H__
/* canonical label */

#include "dg.h"


/* set parameters for the program `nauty' (No AUTomorphism, Yes?) */
#define WORDSIZE        DG_WORDBITS
#define MAXN            DG_NMAX
#define ONE_WORD_SETS 1 /* use one-word set when possible */
#include "nau0s.h"



/* canonical label */
static dg_t *dg_canlabel(dg_t *gout, const dg_t *gin)
{
  static dgword_t g0[DG_NMAX*DG_NW], gc[DG_NMAX*DG_NW];
  static int lab[DG_NMAX], ptn[DG_NMAX], orbits[DG_NMAX];
  static DEFAULTOPTIONS_GRAPH(options);
  static statsblk stats;
#pragma omp threadprivate(g0, gc, lab, ptn, orbits, options, stats)
  int i;
  DGVS_DEFIQ_(iq)
  DG_DEFN_(gin)

  /* nauty uses the highest bit for the first index */
  for (i = 0; i < DG_N_; i++) {
#if DGVS_ONEWORD
    g0[i] = bitrev(gin->c[i]);
#else /* !DGVS_ONEWORD */
    for (iq = 0; iq < DG_NW; iq++)
      g0[i * DG_NW + iq] = bitrev(gin->c[i][iq]);
#endif /* DGVS_ONEWORD */
  }
  options.getcanon = TRUE;
  densenauty(g0, lab, ptn, orbits, &options, &stats, DG_NW, DG_N_, gc);

  gout->n = DG_N_;
  for (i = 0; i < DG_N_; i++) {
#if DGVS_ONEWORD
    gout->c[i] = bitrev(gc[i]);
#else /* !DGVS_ONEWORD */
    for (iq = 0; iq < DG_NW; iq++)
      gout->c[i][iq] = bitrev(gc[i * DG_NW + iq]);
#endif /* DGVS_ONEWORD */
  }
  return gout;
}



typedef struct {
  int n; /* number of vertices */
  int nc; /* number of sets in the partition */
  dgvs_t cs[DG_NMAX]; /* cs[i] is the ith vertex set */
  int cnt[DG_NMAX]; /* cnt[i] is the number of vertices in set i */
} dgpart_t; /* an ordered partition, a sequence of vertex sets */



INLINE void dgpart_unit(dgpart_t *part, int n)
{
  part->n = n;
  part->nc = 1;
  DGVS_MKBITSMASK(part->cs[0], n)
  part->cnt[0] = n;
}



#define dgpart_print(p) dgpart_fprint(p, stdout)
#define dgpart_errprint(p) dgpart_fprint(p, stderr)

/* print out the partition */
static void dgpart_fprint(const dgpart_t *part, FILE *fp)
{
  int ic, k;
  dgvs_t vs;
  dgword_t b;
  DGVS_DEFIQ_(iq)

  for (ic = 0; ic < part->nc; ic++) {
    DGVS_CPY(vs, part->cs[ic])
    while ( dgvs_nonzero(vs) ) {
      DGVS_FIRSTLOW(k, vs, b, iq)
      DGVS_XOR1(vs, b, iq)
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
  dgvs_t all, vs, vs1;

  DGVS_CLEAR(all)
  for (ic = 0; ic < part->nc; ic++) {
    DGVS_CPY(vs, part->cs[ic]);
    /* `all' collects all encountered vertices so far,
     * so if this partition `vs' overlaps with `all'
     * there is a conflict */
    if ( dgvs_nonzero( dgvs_and2(vs1, vs, all) ) ) {
      fprintf(stderr, "cell %d has a previous element %d\n",
          ic, dgvs_first(vs1));
      dgpart_errprint(part);
      return 1;
    }
    sz = dgvs_count(vs); /* count vertices in `vs' */
    DGVS_OR(all, vs) /* add vertices in `vs' to `all' */
    if (sz != part->cnt[ic] || sz == 0) {
      fprintf(stderr, "bad cell %d size %d, or it mismatches %d\n",
          ic, sz, part->cnt[ic]);
      dgpart_errprint(part);
      return 2;
    }
  }
  DGVS_MKBITSMASK(vs, part->n)
  /* `all' must contain all vertices now */
  if ( dgvs_neq(all, vs) ) {
    fprintf(stderr, "partition does not cover all vertices\n");
    dgpart_errprint(part);
    return 3;
  }

  if (equitable) { /* see if the partition is equitable w.r.t. all cells */
    /* loop over the shatterers */
    for (ic = 0; ic < part->nc; ic++) {
      int jc, deg0, deg, k, k0;
      dgvs_t vsi, vsj;
      dgword_t b;
      DGVS_DEFIQ_(iq)

      DGVS_CPY(vsi, part->cs[ic]);
      /* loop over the shatterees */
      for (jc = 0; jc < part->nc; jc++) {
        /* if the cell has only one vertex, it cannot be shattered */
        if (part->cnt[jc] == 1) continue;
        deg0 = k0 = -1;
        /* loop over vertices in the cell jc */
        DGVS_CPY(vsj, part->cs[jc])
        while ( dgvs_nonzero(vsj) ) {
          DGVS_FIRSTLOW(k, vsj, b, iq) /* first vertex `k' in `vsj' */
          DGVS_XOR1(vsj, b, iq) /* remove `b' from `vsj' */
          deg = dg_degvs(g, k, vsi); /* degree of `k' w.r.t. cell `i' */
          if (deg0 < 0) {
            deg0 = deg;
            k0 = k;
          } else if (deg != deg0) {
            fprintf(stderr, "cell %d (point %d) is not equitable w.r.t. cell %d, degs %d vs %d (point %d)\n",
                jc, k, ic, deg, deg0, k0);
            dgpart_print(part);
            dg_print(g);
            return 9;
          }
        }
      }
    }
  }
  return 0;
}



/* return an equitable partition of the graph that is compatible with `g'
 * An equitable partition is a sequence of vertex subsets or cells {Vi},
 * such that for any Vi and Vj (i can be equal to j) and
 * for any vk and vl in Vi, deg(vk, Vj) == deg(vk, Vi) */
INLINE void dg_equipart(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g)
  int i, i0, imax, j, k, deg, ni;
  dgvs_t vsj, vsi;
  static dgvs_t vswdeg[DG_NMAX]; /* vswdeg[i] is the vertex set with degree i */
#pragma omp threadprivate(vswdeg)
  dgword_t b;
  DGVS_DEFIQ_(iq)

  /* loop over shatterers
   * note: the number of cells `part->nc' will be increased later
   * in the loop, so it should not be fixed at the beginning of the loop */
  for (j = 0; j < part->nc; j++) {
    DGVS_CPY(vsj, part->cs[j])
    //printf("shatterer j %d\n", j); dgpart_print(part); getchar();

    /* loop over shatterees
     * In this process, the shatterer `j' is treated as fixed
     * so once the cell `i' is shattered properly, its subcells
     * cannot be further shattered by `j', this means that we
     * only need to loop to the current number of cells imax
     * even if number of cells part->nc grows latter on.
     * Note, the shatterer `j' may be changed during the process.
     * However, if Vj_old cannot shatter a cell i*, but the
     * new Vj_new (a subset of Vj_old) can shatter this cell i*
     * then the subset Vj_old - Vj_new can shatter this cell i* as well
     * And Vj_old - Vj_new contains all newly created subcells appended
     * at the end of the list.
     * Thus, it is okay to jump to the next j even if j is changed */
    imax = part->nc; /* the current number of cells */
    i0 = -1; /* the first shattered cell */
    for (i = 0; i < imax; i++) {
      dgvs_t nzdegs;
      dgword_t bdeg;
      DGVS_DEFIQ_(degq)

      ni = part->cnt[i]; /* number of vertices */
      //die_if (ni == 0, "set %d is empty\n", i);
      if (ni == 1) /* one-vertex cell, cannot be shattered */
        continue;

      /* use count sort to shatter `vsi'
       * sort cells according to their degrees with respect to `vsj'
       * then each nonempty bucket labeled by the degree gives
       * the set of shatterees of the same degree */
      DGVS_CPY(vsi, part->cs[i]);
      /* `nzdegs' has bit `deg' if there is a vertex in cell `vsi'
       * with a degree `deg' w.r.t. `vsj' */
      DGVS_CLEAR(nzdegs)
      /* loop over vertices in the shatteree `vsi' */
      while ( dgvs_nonzero(vsi) ) {
        DGVS_FIRSTLOW(k, vsi, b, iq) /* first vertex `k' in `vsi' */
        DGVS_XOR1(vsi, b, iq) /* remove `b' from `vsi' */
        deg = dg_degvs(g, k, vsj); /* degree of `k' w.r.t. to `vsj' */
        DGVS_MKBIT(deg, bdeg, degq)
        /* vswdeg[deg] already has something */
        if ( !DGVS_HASBIT(nzdegs, bdeg, degq) ) {
          /* mark degree `deg' has been occupied */
          DGVS_OR1(nzdegs, bdeg, degq)
          DGVS_CLEAR(vswdeg[deg]) /* clear the bucket for `deg' */
        }
        DGVS_XOR1(vswdeg[deg], b, iq) /* add `b' into the bucket `vswdeg[deg]' */
      }
      /* if there is only one occupied cell, no shattering */
      if ( dgvs_count(nzdegs) == 1) continue;
      /* now nonzero vswdeg[deg] collects a subset of `vsi' with degree `deg'
       * we loop over occupied degrees */
      k = 0; /* number of shattered cells */
      while ( dgvs_nonzero(nzdegs) ) {
        DGVS_FIRSTLOW(deg, nzdegs, bdeg, degq)
        DGVS_XOR1(nzdegs, bdeg, degq) /* remove `bdeg' from `nzdegs' */
        //printf("shatter %d k %d deg %d, part->nc %d\n", i, k, deg, part->nc);
        //dgvs_printn(vswdeg[deg], "vs[deg]");
        if (k == 0) { /* replace cell `i' by the first partition */
          DGVS_CPY(part->cs[i], vswdeg[deg])
          part->cnt[i] = dgvs_count( part->cs[i] );
          k = part->nc; /* the next set begins here */
        } else {
          DGVS_CPY(part->cs[k], vswdeg[deg])
          part->cnt[k] = dgvs_count( part->cs[k] );
          part->nc = ++k;
          if (i0 < 0) i0 = i; /* update i0, if cell i has been shattered */
        }
      }
      //dgpart_print(part); getchar();
      /* if the partition has n cells, return */
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



/* check if a permutation of n is good, return 0 on success */
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
 * if `recur' is true, try to artificially break the equitable partition
 *   until it becomes discrete
 * the speed of using recur or not is very similar, recommended */
INLINE void dg_permequipart(int *perm, const dg_t *g, int recur)
{
  DG_DEFN_(g)
  dgpart_t part;
  dgvs_t vs;
  int ic, k, i;
  dgword_t b;
  DGVS_DEFIQ_(iq)

  /* initially set part as the unit partition */
  dgpart_unit(&part, DG_N_);

  if (recur) {
    /* recursively call dg_equipart() until part becomes discrete */
    while (part.nc < DG_N_) {
      /* call dg_equipart() the shatter `part' */
      dg_equipart(&part, g);

      if (part.nc == DG_N_) break;
      /* find the first unfixed cell */
      for (ic = 0; ic < part.nc; ic++)
        if (part.cnt[ic] > 1)
          break;
      /* artificially break the first cell */
      DGVS_CPY(vs, part.cs[ic])
      DGVS_FIRSTBIT(vs, b, iq) /* extract the first bit `b' from `vs' */
      DGVS_XOR1(vs, b, iq) /* remove `b' from `vs' */

      DGVS_CPY(part.cs[part.nc], vs) /* append a cell with the new `vs' */
      part.cnt[part.nc] = part.cnt[ic] - 1;

      DGVS_ONEBIT(part.cs[ic], b, iq) /* change this cell to a single bit */
      part.cnt[ic] = 1;

      part.nc += 1;
    }
    /* now each cell is discrete, and the partition corresponds to a permutation */
    for (ic = 0; ic < DG_N_; ic++) {
      perm[ dgvs_first(part.cs[ic]) ] = ic;
    }
  } else { /* single pass */
    dg_equipart(&part, g);
    /* find the first permutation compatible with the partitions */
    /* loop over cells in the partition */
    for (i = 0, ic = 0; ic < part.nc; ic++) {
      DGVS_CPY(vs, part.cs[ic])
      while ( dgvs_nonzero(vs) ) {
        DGVS_FIRSTLOW(k, vs, b, iq) /* first bit `b' in vs */
        DGVS_XOR1(vs, b, iq) /* remove bit `b' from `vs' */
        perm[ k ] = i++;
      } /* loop over vertices in a partition */
    } /* loop over partitions */
  } /* recur or single-pass */
}



/* get a permutation `perm' compatible with the degree sequence for graph `g' */
INLINE void dg_permdegseq(int *perm, const dg_t *g)
{
  DG_DEFN_(g)
  int i, k, deg;
  dgvs_t cvs[DG_NMAX], vs;
  dgword_t b;
  DGVS_DEFIQ_(iq)
#define DG_DEGMIN 0 /* can be 2 for biconnected graphs */

  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    DGVS_CLEAR( cvs[deg] )
  }
  /* use count sort to collect vertices with the same degrees */
  for (i = 0; i < DG_N_; i++) {
    /* cvs[d] collects vertices (as bits)
     * with the same degree, given by dg_deg(g, i) */
    DGVS_ADD( cvs[dg_deg(g, i)], i )
  }

  i = 0; /* rank of the current vertex, sorted by the degree */
  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    DGVS_CPY(vs, cvs[deg]) /* copy the vertex set */
    /* now `vs' collects all vertices with the same degree, `deg' */
    while ( dgvs_nonzero(vs) ) { /* loop over vertices in this set `vs' */
      DGVS_FIRSTLOW(k, vs, b, iq) /* extract the first bit `b' from `vs' */
      DGVS_XOR1(vs, b, iq) /* remove `b' from `vs' */
      perm[ k ] = i++; /* vertex `k' is the ith lowest-degree vertex */
    }
  }
}



/* get a partition in which each cell contains all vertices of the same degree */
INLINE void dg_partdegseq(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g)
  int i, deg, ic;
  dgvs_t cvs[DG_NMAX];
#define DG_DEGMIN 0 /* can be 2 for biconnected graphs */

  for (i = DG_DEGMIN; i < DG_N_; i++) {
    DGVS_CLEAR( cvs[i] )
  }
  /* use count sort to collect vertices with the same degrees
   * `cvs[deg]' collects vertices of the same degree `deg' */
  for (i = 0; i < DG_N_; i++) {
    DGVS_ADD(cvs[dg_deg(g, i)], i)
  }
  /* construct a partition accordingly */
  part->n = g->n;
  ic = 0;
  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    dgvsref_t vs = cvs[deg];
    if ( !dgvs_nonzero(vs) ) continue;
    DGVS_CPY(part->cs[ic], vs)
    part->cnt[ic] = dgvs_count(vs);
    ic++;
  }
  part->nc = ic;
}



/* rearrange vertices in `gin' according to the permutation `perm'
 * save the resulting graph in `gout' */
INLINE dg_t *dg_permute(dg_t * RESTRICT gout, const dg_t * RESTRICT gin, const int *perm)
{
  DG_DEFN_(gin)
  int i, j;
  dgvs_t ci, ci_out;

  gout->n = DG_N_;
  dg_empty(gout); /* adding this appears to make the function faster? */
  for (i = 0; i < DG_N_; i++) {
    DGVS_CPY(ci, gin->c[i])
    DGVS_CLEAR(ci_out)
    for (j = 0; j < DG_N_; j++) {
      if ( DGVS_HAS(ci, j) ) {
        DGVS_ADD(ci_out, perm[j])
      }
    }
    DGVS_CPY(gout->c[ perm[i] ], ci_out);
  }
  return gout;
}



/* check if a is a vertex permutation of b */
INLINE int dg_checkiso(const dg_t *a, const dg_t *b, const int *perm)
{
  DG_DEFN_(b)
  int i, j;

  for (i = 0; i < DG_N_; i++) {
    for (j = 0; j < DG_N_; j++) {
      if (DGVS_HAS(a->c[i], j)
       != DGVS_HAS(b->c[perm[i]], perm[j]))
        return 1;
    }
  }
  return 0;
}



/* get a representative isomorphic graph `gout' of the input graph `gin'
 * `level' represent the depth of the search
 * if level == 4, it returns the unique canonical label of the graph
 *   all isomorphic graphs share the same canonical label
 * if level == 3, it returns the first node from McKay's search of the
 *    canonical label
 * if level == 2, it returns a graph that is compatible with equitable
 *    partition (hierarchical degree sequence) of the graph
 * if level == 1, it returns a graph compatible with the degree sequence
 *
 * In all cases, `gout' is isomorphic to `gin' since the former is produced
 *  from some vertex permutation of the latter.
 * But isomorphic graphs may not be mapped to the same `gout' if level < 3
 * with a higher `level' the redundancy reduces, but the speed goes down.
 * gout->n will be proper set here
 * */
INLINE dg_t *dg_repiso(dg_t *gout, const dg_t *gin, int level)
{
  int perm[DG_NMAX + 1];

  if (level == 0) {
    /* disable isomorphism, simply copy the graph */
    DG_DEFN_(gin)
    gout->n = DG_N_;
    return dg_copy(gout, gin);
  }

  if (level >= 4 || level < 0) { /* return the canonical label */
    /* use NAUTY to find the unique isomorphism */
    return dg_canlabel(gout, gin);
  }

  /* compute a permutation `perm' of the vertices according to the `level'
   * the permutation represents an isomorphic graph to `g'
   * with a higher `level', there are fewer possible permutations or isomorphisms
   * but it takes longer to compute it */
  if (level == 1) { /* degree sequence */
    dg_permdegseq(perm, gin);
  } else if (level == 2) { /* equitable partition */
    dg_permequipart(perm, gin, 0);
  } else if (level == 3) { /* equitable partition */
    dg_permequipart(perm, gin, 1);
  }
  //die_if (dg_checkperm(perm, gin->n) != 0, "bad permutation level %d\n", level);

  /* compute the isomorphic graph `gout' of the input graph `gin'
   * according to the permutation `perm' computed above */
  dg_permute(gout, gin, perm);
  //die_if (dg_checkiso(gout, gin, perm) != 0, "bad isomorphism!\n");
  return gout;
}



/* enumerate all permutations compatible with a partition  */
static int dg_part2permls(int *perms, int nmax, dgpart_t *part)
{
  int i, ic, k, pcnt, n = part->n, nc = part->nc;
  dgvs_t vs, vs1, unused;
  dgword_t b;
  int st[DG_NMAX], top = 0;
  //int stdefid[DG_NMAX];
  dgvs_t stvs[DG_NMAX], stavail[DG_NMAX];
  DGVS_DEFIQ_(iq)

  /* vertices available in each stack position */
  for (i = 0, ic = 0; ic < nc; ic++) {
    //int defid = (part->cnt[ic] == 1) ? dgvs_first(part->cs[ic]) : -1;
    for (k = 0; k < part->cnt[ ic ]; k++, i++) {
      /* vertices in the cell `ic' share the same vertex set */
      DGVS_CPY(stvs[i], part->cs[ ic ]);
      //stdefid[i] = defid;
    }
  }
  //dgpart_print(part);
  //for (i = 0; i < n; i++) printf("%d %d\n", i, stdefid[i]);
  //die_if (i != n, "bad partition i %d != n %d\n", i, n);

  pcnt = 0;
  top = 0;
  DGVS_MKBITSMASK(unused, n)
  DGVS_CPY(stavail[top], stvs[top])
  while (1) {
    /* determine the vertex on the top */
    DGVS_CPY(vs, stavail[top])
    if ( vs ) { /* we still have vertices on this level */
      DGVS_FIRSTLOW(k, vs, b, iq) /* select the vertex `k' from `vs' */
      st[top] = k; /* set the vertex on the top as `k' */
      //printf("push top %d, k %d, stavail %s, unused %s\n", top, k, dg_strset1(stavail[top]), dg_strset1(unused) );
      /* remove `b' from the available set `stavail[top]'
       * equivalent to: stavail[top] = vs ^ b */
      DGVS_XOR1(vs, b, iq)
      DGVS_CPY(stavail[top], vs)
      if (top < n - 2) { /* push */
        dgvsref_t vsr;
        DGVS_XOR1(unused, b, iq) /* `unused' is the current vertex set
                                    with the current stack pointer
                                    remove `b' from this set */
        /* if the next level `top + 1' share the same cell with the
         * current level, we update the available vertex set
         * by excluding used vertices */
        vsr = stvs[top + 1]; /* all vertices in the cell */
        if ( dgvs_eq(vsr, stvs[top]) ) { /* same cell */
          DGVS_AND2(stavail[top + 1], vsr, unused) /* copy the unused list */
        } else { /* different cells, copy all vertices in the cell */
          DGVS_CPY(stavail[top + 1], vsr)
        }
        top++;
        continue;
      }
      /* here top == n - 2, this means that we are at the next to
       * the top level, then no need to go up further */
      /* determine the last vertex, which is given by the
       * single bit in the vertex set given by vs1 = unused ^ b */
      DGVS_CPY(vs1, unused)
      DGVS_XOR1(vs1, b, iq)
      //die_if (dgvs_iszero(vs1), "no vertex on the last level, top %d\n", top);
      st[top + 1] = dgvs_bit2id(vs1);
      //dgpart_print(part); dg_printperm(st, n); getchar();
      /* TODO: remove this */
      //die_if (dg_checkperm(st, n) != 0, "bad permutation n %d\n", n);
      /* save the permutation */
      for (k = 0; k < n; k++)
        perms[pcnt * n + st[k]] = k;
      if (++pcnt >= nmax) return pcnt;
      /* do not increase top, just explore other choices on this level */
    } else { /* no vertex left on this level, pop */
      if (--top < 0) break;
      DGVS_ADD(unused, st[top]);
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
  dgpart_t part = {0};
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
  /* obtain a list of possible permutations that are compatible
   * with the vertex partition `part' */
  ncnt = dg_part2permls(dgaut_perms_, nmax, &part);
  /* convert each permutation to a code */
  for (ic = 0, i = 0; i < ncnt; i++) {
    /* construct the graph `ng' from `g' by permutation `i' */
    dg_permute(ng, g, dgaut_perms_ + i * n);
    /* encode the graph `ng' into `codes[ic]' */
    dg_encode(ng, codes + ic * cwords);
    /* avoid duplicated codes see if the code is new */
    for (j = 0; j < ic; j++)
      if (DG_CEQ(codes + ic * cwords, codes + j * cwords, cwords))
        break;
    /* increment ic only if the code is new */
    if (j >= ic) ic++;
  }
  return ic;
}



#endif /* DGAUT_H__ */

