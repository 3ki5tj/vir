#ifndef DGMAPL_H__
#define DGMAPL_H__

#ifndef DGMAPL_NMAX
#define DGMAPL_NMAX 9
#endif

#if DGMAPL_NMAX >= 10 /* for future extension, if any */
#ifndef CODEBITS
#define CODEBITS 64 /* for the table key is over 2^32 */
#endif /* CODEBITS */
#endif /* DGMAPL_NMAX >= 10 */

#include "dg.h"
#include "dgsc.h"
#if DGMAPL_NMAX <= 9
#ifndef RJW32
#define RJW32 1 /* 32-bit integer is good enough for the
                   RJW recursion at n = 9 < 14, and it is faster
                   on a 32-bit machine */
#endif /* RJW32 */
#endif /* DGMAPL_NMAX <= 9 */
#include "dgrjw.h"
#include "dgring.h"


#ifdef DGMAPL_INT32
/* if DGMAPL_NMAX >= 10, 32-bit integer may be preferrable */
typedef int32_t dgmapl_int_t; /* data type for fb and nr */
#define DGMAPL_MAX ((dgmapl_int_t)  2147483647) /* maximal allowable 32-bit integer */
#define DGMAPL_BAD ((dgmapl_int_t) -2147483647) /* minimal 32-bit integer */
#else
typedef int16_t dgmapl_int_t; /* data type for fb and nr */
#define DGMAPL_MAX ((dgmapl_int_t)  32767) /* maximal allowable 16-bit integer */
#define DGMAPL_BAD ((dgmapl_int_t) -32767) /* minimal 16-bit integer */
#endif



/* encode a graph according to the permuation `st' */
INLINE code_t dgmapl_encode(const dg_t *g, int k, int *st)
{
  int n = g->n, i, j, jmax, pid;
  code_t c;

  for (pid = 0, c = 0, i = 2; i < n; i++) {
    jmax = (i <= k) ? (i - 1) : i;
    for (j = 0; j < jmax; j++, pid++)
      if ( dg_linked(g, st[i], st[j]) )
        c |= (code_t) 1 << pid;
  }
  return c;
}



/* find a chain of k links
 * set `ring' to 1 to find a closed ring */
INLINE code_t dgmapl_getchain(const dg_t *g, int k, int *st)
{
  int n = g->n, err = 1, top, i;
  code_t vs, c, b, ms[DGMAPL_NMAX + 1];

  top = 0;
  ms[top] = vs = mkbitsmask(n);
  while (1) {
    if (ms[top] != 0) { /* push */
      BITFIRSTLOW(st[top], ms[top], b);
      vs ^= b;
      if (top >= k) {
        err = 0;
        goto LOOPEND;
      }
      /* remaining vertices and neighbors of top */
      ms[top + 1] = vs & g->c[ st[top] ];
      top++;
    } else { /* pop */
      if (--top < 0) break;
      b = MKBIT(st[top]);
      vs ^= b; /* add b back to the unused vertex list */
      ms[top] ^= b; /* remove b from ms[top] */
    }
  }

LOOPEND:
  if (err) return 0;
  /* settle the indices of the rest vertices */
  for (i = k + 1, c = vs; c; c ^= b, i++)
    BITFIRSTLOW(st[i], c, b);
  die_if (i != n, "i %d, n %d, vs %#x, k %d\n", i, n, (unsigned) vs, k);
  return dgmapl_encode(g, k, st);
}



/* find all chains of k links in the graph g, and assign the values
 * of the graphs from permutating indices to val
 * `st' will be destoryed */
INLINE int dgmapl_save2full(dgmapl_int_t (*arr)[2], const dgmapl_int_t val[2],
    const dg_t *g, int k, int top, int *st)
{
#define DGMAPL_LISTSIZE 200
  int n = g->n, cnt = 0, i;
  code_t vs, c, b;
  code_t ms[DGMAPL_NMAX + 1];
  code_t ls[DGMAPL_LISTSIZE]; /* buffer to save updated values */
  int lscnt = 0, j;

  vs = ms[0] = mkbitsmask(n);
  if (top != 0) {
    die_if (k != top, "top %d must be %d\n", top, k);
    /* reconstruct the stack from st[] */
    for (i = 0; i < top; i++) {
      b = MKBIT(st[i]);
      vs ^= b; /* remove b from the stack */
      ms[i + 1] = vs & g->c[ st[i] ];
      if ( !dg_linked(g, st[i], st[i+1]) ) {
        fprintf(stderr, "i %d, no edge between %d and %d\n", i, st[i], st[i+1]);
        for (i = 0; i <= top; i++) printf("%d, %d\n", i, st[i]);
        dg_print(g);
        exit(1);
      }
    }
  }
  //dg_print(g);
  //for (ms[top+1] = 0, i = 0; i <= top; i++) printf("%d %#5x %#5x\n", st[i], (unsigned) ms[i], (unsigned) MKBIT(st[i]));
  //printf("initial top %d, vs %#5x\n\n", top, (unsigned) vs);
  while (1) {
    if (ms[top] != 0) BITFIRSTLOW(st[top], ms[top], b);
    //printf("top %d, ms %#5x, st %d, b %#5x, vs %#5x\n", top, (unsigned) ms[top], st[top], (unsigned) MKBIT(st[top]), (unsigned) vs); getchar();
    if (ms[top] != 0) { /* push */
      BITFIRSTLOW(st[top], ms[top], b);
      ms[top + 1] = vs ^ b;
      if (top < k) /* enforcing the connectivity top-(top+1) */
        ms[top + 1] &= g->c[ st[top] ];
      if (top + 1 < n) { /* to the next step to avoid popping */
        vs ^= b;
        top++;
        continue;
      }
      c = dgmapl_encode(g, k, st);
      //for (i = 0; i < n; i++) printf("%d %#5x %#5x\n", st[i], (unsigned) ms[i], (unsigned) MKBIT(st[i]) );
      //printf("above is %d, code %#8x\n\n", cnt, (unsigned) c); //getchar();
      if ( (arr[c][0] != DGMAPL_BAD && arr[c][0] != val[0])
        || (arr[c][1] != DGMAPL_BAD && arr[c][1] != val[1]) ) {
        fprintf(stderr, "%#" PRIx64 " has been occupied, %d,%d (new) vs %d,%d (old)\n",
            (uint64_t) c, val[0], val[1], arr[c][0], arr[c][1]);
        dg_print(g);
        exit(1);
      }

      if (arr[c][0] == DGMAPL_BAD && arr[c][1] == DGMAPL_BAD) {
        ls[lscnt] = c;
        if (++lscnt >= DGMAPL_LISTSIZE) {
#pragma omp critical
          { /* dump */
            for (j = 0; j < lscnt; j++) {
              arr[ls[j]][0] = val[0];
              arr[ls[j]][1] = val[1];
            }
          }
          lscnt = 0; /* clear the buffer */
        }
      }
      cnt++;
      /* fall through to pop */
    }
    /* pop */
    if (--top < 0) break;
    b = MKBIT(st[top]);
    vs ^= b;
    ms[top] ^= b; /* remove b from ms[top] */
  }
  //printf("got %d\n", cnt); getchar();
#pragma omp critical
  { /* dump */
    for (j = 0; j < lscnt; j++) {
      arr[ls[j]][0] = val[0];
      arr[ls[j]][1] = val[1];
    }
  }
  return cnt;
}



/* for n = 9, nr <= 9!/18 = 20160, |fb| <= 7! = 5040
 * so they can be contained in a 16-bit integer */
typedef struct {
  int n;
  int k;
  dgmapl_int_t (*fbnr)[2];
  dgmapl_int_t *nr;
} dgmapl_t;



INLINE dgmapl_t *dgmapl_open(int n, int k)
{
  dgmapl_t *mapl;
  uint64_t i, size;
  static const int kdef[] = {0,
    0, 1, 2, 3, 4, 4, 5, 6, 8, 9};

  die_if (n > DGMAPL_NMAX, "n %d is too large\n", n);
  xnew(mapl, 1);
  mapl->n = n;
  if (k <= 0) /* choose a default k */
    k = kdef[n];
  mapl->k = k;
  size = (uint64_t) 1u << (n * (n - 1) / 2 - k);
  xnew(mapl->fbnr, size);
  for (i = 0; i < size; i++) {
    mapl->fbnr[i][0] = DGMAPL_BAD;
    mapl->fbnr[i][1] = DGMAPL_BAD;
  }
  return mapl;
}



INLINE void dgmapl_close(dgmapl_t *mapl)
{
  free(mapl->fbnr);
  free(mapl);
}



#define dg_fbnr_Lookup(g, k, nr) dg_fbnr_Lookup0(g, k, nr, 0, NULL, NULL)

/* compute fb and nr by the larger lookup table
 * set k = 0 to use the default search length */
INLINE double dg_fbnr_Lookup0(const dg_t *g, int k, double *nr,
    int nocsep, int *ned, int *degs)
{
#ifdef DGMAPL_DEBUG
  static uint64_t mis = 0, tot = 0, cnt = 0, hit = 0;
#endif
  static dgmapl_t *mapl[DGMAPL_NMAX + 1];
  static int st[DGMAPL_NMAX + 1];
#pragma omp threadprivate(st)

  code_t c;
  int n = g->n, hasnew;
  dgmapl_int_t ifbnr[2];
  dgmapl_t *tb;
  double fb;

  die_if (n <= 0 || n > DGMAPL_NMAX, "bad n %d\n", n);

  /* allocate memory */
  if (mapl[n] == NULL) {
#ifdef _OPENMP
    /* We test mapl[n] twice here. If mapl[n] != NULL, then the map
     * memory has been allocated, and we may safely skip the critical
     * block.  If mapl[n] == NULL, the if clause inside the critical
     * block ensures that the memory is not allocated twice.  That is,
     * after the first thread leaves the critical block, the second
     * thread can see `mapl[n]' is no longer NULL, and skip the block */
#pragma omp critical
    {
      if (mapl[n] == NULL) {
#endif /* _OPENMP */
        clock_t t0 = clock();
        mapl[n] = dgmapl_open(n, k);
        fprintf(stderr, "%4d: initialized fb/nr-mapl for %d, k %d, time %gs\n",
            inode, n, mapl[n]->k, 1.*(clock() - t0)/CLOCKS_PER_SEC);
#ifdef _OPENMP
      }
    } /* omp critical */
#endif /* _OPENMP */
  }

  tb = mapl[n];
  c = dgmapl_getchain(g, tb->k, st);
  if (c == 0) { /* the lookup table is not applicable, do it directly */
#ifdef DGMAPL_DEBUG
#pragma omp critical
    { mis += 1; tot += 1; }
#endif
    *nr = 0;
    fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
    return fb;
  }

  ifbnr[0] = tb->fbnr[c][0];
  ifbnr[1] = tb->fbnr[c][1];
  hasnew = 0; /* if we have computed new values */
  if (ifbnr[0] != DGMAPL_BAD) {
#ifdef DGMAPL_DEBUG
#pragma omp atomic
    hit += 1;
#endif
    fb = ifbnr[0];
  } else {
    fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
    if (fabs(fb) < DGMAPL_MAX) {
      ifbnr[0] = (dgmapl_int_t) ((fb < 0) ? (fb - .5) : (fb + .5));
      hasnew = 1;
    }
  }

  if (ifbnr[1] != DGMAPL_BAD) {
    *nr = ifbnr[1];
  } else {
    *nr = dg_nring_mixed0(g, ned, degs);
    if (*nr < DGMAPL_MAX) {
      ifbnr[1] = (dgmapl_int_t) (*nr + .5);
      hasnew = 1;
    }
  }

  if (hasnew) {
#ifdef DGMAPL_SIMPLEUPDATE
#pragma omp critical
    {
      tb->fbnr[c][0] = ifbnr[0];
      tb->fbnr[c][1] = ifbnr[1];
#ifdef DGMAPL_DEBUG
      cnt += 1;
#endif
    }
#else /* full update, slower but more effective */
#ifdef DGMAPL_DEBUG
    cnt +=
#endif
      /* save fb and nr for different ways of tracing the graph */
      dgmapl_save2full(tb->fbnr, ifbnr, g, tb->k, tb->k, st);
#endif /* DGMAPL_SIMPLEUPDATE */
  }

#ifdef DGMAPL_DEBUG /* debugging code */
#pragma omp critical
  {
    tot += 1;
#ifndef DGMAPL_NSTREP
#define DGMAPL_NSTREP 1000000
#endif
    if (tot % DGMAPL_NSTREP == 0)
      fprintf(stderr, "%4d: cnt %g, hits %g/%g = %g%%, misses %g/%g = %g%% graphs\n",
          inode, 1.*cnt, 1.*hit, 1.*tot, 100.*hit/tot, 1.*mis, 1.*tot, 100.*mis/tot);
  }
#endif /* DGMAPL_DEBUG */

  return fb;
}




#endif

