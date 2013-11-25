/* get a permutation compatible with the degree sequence */
INLINE void dg_permdegseq0(int *perm, const dg_t *g)
{
  DG_DEFN_(g);
  int i, k, deg;
  code_t cvs[DG_NMAX], vs, b;
#define DG_DEGMIN 0 /* can be 2 for biconnected graphs */

  for (i = DG_DEGMIN; i < DG_N_; i++) cvs[i] = 0;
  /* use count sort to collect vertices with the same degrees */
  for (b = 1, i = 0; i < DG_N_; i++, b <<= 1)
    cvs[dg_deg(g, i)] |= b;
  i = 0;
  for (deg = DG_DEGMIN; deg < DG_N_; deg++) {
    for (vs = cvs[deg]; vs; vs ^= b) { /* loop over this partition */
      BITFIRSTLOW(k, vs, b);
      perm[k] = i++;
    }
  }
}



/* get a permutation compatible with the degree sequence */
INLINE void dg_permdegseq1(int *perm, const dg_t *g)
{
  DG_DEFN_(g);
  int i, k, deg;
  code_t cvs[DG_NMAX], vs, b;
  code_t nzdegs, nz, bdeg;

  /* use count sort to collect vertices with the same degrees */
  nzdegs = 0;
  for (i = 0; i < DG_N_; i++) {
    deg = dg_deg(g, i);
    b = MKBIT(i);
    bdeg = MKBIT(deg);
    if (bdeg & nzdegs) { /* cvs[deg] has been initialized */
      cvs[deg] |= b;
    } else {
      cvs[deg] = b;
      nzdegs |= bdeg; /* mark degree `deg' has been occupied */
    }
  }
  i = 0;
  /* loop over occupied degrees */
  for (nz = nzdegs; nz; nz ^= bdeg) {
    BITFIRSTLOW(deg, nz, bdeg);
    for (vs = cvs[deg]; vs; vs ^= b) { /* loop over this partition */
      BITFIRSTLOW(k, vs, b);
      die_if (i > DG_N_, "i %d is too large\n", i);
      perm[k] = i++;
    }
  }
}




/* return an equitable partition of the graph that is compatible with part
 * An equitable partition is a sequence of vertex subsets or cells {Vi},
 * such that for any Vi and Vj (i can be equal to j) and
 * for any vk and vl in Vi, deg(vk, Vj) == deg(vk, Vi) */
INLINE void dg_equipart(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g);
  int i, i0, imax, j, k, deg, ni;
  code_t vsj, vsi, vs, b;
  code_t vswdeg[DG_NMAX]; /* vswdeg[i] is the vertex set with degree i */

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
      die_if (ni == 0, "set %d is empty\n", i);
      if (ni == 1) { /* one-vertex cell, nothing to shatter */
        continue;
      } else if (ni == 2) { /* cell has two vertices */
        int k2, deg2;
        code_t b2;

        BITFIRSTLOW(k, vsi, b); /* get the first vertex */
        deg = dg_degvs(g, k, vsj);
        b2 = vsi ^ b; /* get the second vertex */
        k2 = BIT2ID(b2);
        deg2 = dg_degvs(g, k2, vsj);
        if (deg == deg2) { /* degrees are the same, cannot shatter */
          continue;
        } else if (deg < deg2) {
          part->cs[i] = b;
          part->cnt[i] = 1;
          part->cs[part->nc] = b2;
          part->cnt[part->nc] = 1;
        } else {
          part->cs[i] = b2;
          part->cnt[i] = 1;
          part->cs[part->nc] = b;
          part->cnt[part->nc] = 1;
        }
        if (i0 < 0) i0 = i;
        part->nc += 1;
      } else {
        /* general case, we use count sort to shatter vsi, that is,
         * we sort cells according to their degrees with respect to vsj
         * then each nonempty box labeled by the degree gives the shatterer */
        for (k = 0; k < DG_N_; k++) vswdeg[k] = 0;
        /* loop over vertices in the shatteree */
        for (vs = vsi; vs; vs ^= b) {
          BITFIRSTLOW(k, vs, b);
          deg = dg_degvs(g, k, vsj);
          vswdeg[deg] |= b;
        }
        /* now nonzero vswdeg[x] collects a subset of vsi with deg x */
        for (k = 0, deg = 0; deg < DG_N_; deg++) {
          if (vswdeg[deg] == 0) continue;
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




/* return an equitable partition of the graph that is compatible with part
 * An equitable partition is a sequence of vertex subsets or cells {Vi},
 * such that for any Vi and Vj (i can be equal to j) and
 * for any vk and vl in Vi, deg(vk, Vj) == deg(vk, Vi) */
INLINE void dg_equipart1(dgpart_t *part, const dg_t *g)
{
  DG_DEFN_(g);
  int i, i0, imax, j, k, deg, ni;
  code_t vsj, vsi, vs, b;
  code_t vswdeg[DG_NMAX]; /* vswdeg[i] is the vertex set with degree i */

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
      die_if (ni == 0, "set %d is empty\n", i);
      if (ni == 1) { /* one-vertex cell, nothing to shatter */
        continue;
      } else {
        code_t nzdegs, nz, bdeg;
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







