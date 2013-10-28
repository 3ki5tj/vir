#ifndef DGMAPL_H__
#define DGMAPL_H__

#ifndef DGMAPL_NMAX
#define DGMAPL_NMAX 9
#endif

#if DGMAPL_NMAX >= 10
#ifndef CODEBITS
#define CODEBITS 64 /* for the table key is over 2^32 */
#endif /* CODEBITS */
#endif /* DGMAPL_NMAX >= 10 */

#if defined(N) && N > DGMAPL_NMAX
  #define DGMAPL_EXISTS 0
#else
  #define DGMAPL_EXISTS 1
#endif


#include "dg.h"
#include "dgsc.h"
#include "dgrjw.h"
#include "dgring.h"



typedef int32_t dgmapl_int_t; /* integer type to contain fb and nr upto DGMAPL_NMAX */

/* a short integer of three bytes */
typedef unsigned char byte3_t[3];


#ifdef DGMAPL_INT32
/* if DGMAPL_NMAX >= 10, 32-bit integer may be preferrable */
typedef int32_t dgmapl_fb_t; /* internal type to save fb and nr */
#define DGMAPL_MAX ((dgmapl_int_t) 0x7f7f7f7f)
#define DGMAPL_BAD ((dgmapl_int_t) 0x80808080)

#elif defined(DGMAPL_INT24) /* 3-byte, or 24-bit, integer */
/* this is only needed for DGMAPL_NMAX = 10, and in this case,
 * it is only effective for D == 2, and possibly 3, but not higher */
typedef byte3_t dgmapl_fb_t;
#define DGMAPL_MAX ((dgmapl_int_t) 0x7f7f7f)
#define DGMAPL_BAD ((dgmapl_int_t) 0xff808080)

#else /* 2-byte, or 16-bit, integer, default for DGMAPL_NMAX < 64 */
typedef int16_t dgmapl_fb_t;
#define DGMAPL_MAX ((dgmapl_int_t) 0x7f7f)
#define DGMAPL_BAD ((dgmapl_int_t) 0xffff8080)
#endif



INLINE unsigned char *itob3(uint32_t u, byte3_t b)
{
  b[0] = (unsigned char) (u & 0xff);
  b[1] = (unsigned char) ((u >> 8) & 0xff);
  b[2] = (unsigned char) ((u >> 16) & 0xff);
  return b;
}



INLINE int b3toi(const byte3_t b)
{
  unsigned i = b[0] | ((unsigned) b[1] << 8) | ((unsigned) b[2] << 16);
  if (b[2] & 0x80) i |= ((unsigned) 0xff << 24);
  return (int) i;
}



#ifdef DGMAPL_INT24
#define DGMAPL_GETI(x)      b3toi(x)
#define DGMAPL_SETI(x, i)   itob3(i, x)

#else /* 16-bit or 32-bit integer */
#define DGMAPL_GETI(x)      (x)
#define DGMAPL_SETI(x, i)   ((x) = (dgmapl_fb_t) (i))

#endif



INLINE void dgmapl_geti2(dgmapl_int_t * RESTRICT i, dgmapl_fb_t * RESTRICT x)
{
  i[0] = DGMAPL_GETI(x[0]);
  i[1] = DGMAPL_GETI(x[1]);
}



INLINE void dgmapl_seti2(dgmapl_fb_t * RESTRICT x, dgmapl_int_t * RESTRICT i)
{
  DGMAPL_SETI(x[0], i[0]);
  DGMAPL_SETI(x[1], i[1]);
}



/* encode a graph according to the permuation `st' */
INLINE code_t dgmapl_encode(const dg_t *g, int k, int *st)
{
  int i, j, jmax, pid;
  DG_DEFN_(g);
  code_t c;

  for (pid = 0, c = 0, i = 2; i < DG_N_; i++) {
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
  int err = 1, top, i;
  DG_DEFN_(g);
  DG_DEFMASKN_();
  code_t vs, c, b, ms[DGMAPL_NMAX + 1];

  top = 0;
  ms[top] = vs = DG_MASKN_;
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
  die_if (i != DG_N_, "i %d, n %d, vs %#x, k %d\n", i, DG_N_, (unsigned) vs, k);
  return dgmapl_encode(g, k, st);
}



/* find all chains of k links in the graph g, and assign the values
 * of the graphs from permutating indices to val
 * `st' will be destoryed */
INLINE int dgmapl_save2full(dgmapl_fb_t (* RESTRICT arr)[2],
    dgmapl_int_t * RESTRICT val, const dg_t * RESTRICT g,
    int k, int top, int * RESTRICT st)
{
#define DGMAPL_LISTSIZE 200
  int cnt = 0, i;
  DG_DEFN_(g);
  DG_DEFMASKN_();
  code_t vs, c, b;
  code_t ms[DGMAPL_NMAX + 1];
  code_t ls[DGMAPL_LISTSIZE]; /* buffer to save updated values */
  int lscnt = 0, j;
  dgmapl_int_t iarr[2];

  vs = ms[0] = DG_MASKN_;
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
    //if (ms[top] != 0) BITFIRSTLOW(st[top], ms[top], b);
    //printf("top %d, ms %#5x, st %d, b %#5x, vs %#5x\n", top, (unsigned) ms[top], st[top], (unsigned) MKBIT(st[top]), (unsigned) vs); getchar();
    if (ms[top] != 0) { /* push */
      BITFIRSTLOW(st[top], ms[top], b);
      ms[top + 1] = vs ^ b;
      if (top < k) /* enforcing the connectivity top-(top+1) */
        ms[top + 1] &= g->c[ st[top] ];
      if (top + 1 < DG_N_) { /* to the next step to avoid popping */
        vs ^= b;
        top++;
        continue;
      }
      c = dgmapl_encode(g, k, st);
      //for (i = 0; i < n; i++) printf("%d %#5x %#5x\n", st[i], (unsigned) ms[i], (unsigned) MKBIT(st[i]) );
      //printf("above is %d, code %#8x\n\n", cnt, (unsigned) c); //getchar();
      dgmapl_geti2(iarr, arr[c]);
      if ( (iarr[0] != DGMAPL_BAD && iarr[0] != val[0])
        || (iarr[1] != DGMAPL_BAD && iarr[1] != val[1]) ) {
        fprintf(stderr, "%#" PRIx64 " has been occupied, %d,%d (new) vs %d,%d (old)\n",
            (uint64_t) c, val[0], val[1], iarr[0], iarr[c]);
        dg_print(g);
        exit(1);
      }

      /* we update the table only if both fb and nr are bad */
      if ( iarr[0] == DGMAPL_BAD && iarr[1] == DGMAPL_BAD ) {
#ifdef DGMAPL_NODUP /* avoid duplicates for slow memory */
        /* see if we have already have c */
        for (j = 0; j < lscnt; j++) if (ls[j] == c) break;
#else
        j = lscnt;
#endif
        if (j == lscnt) { /* c is new to the list */
          ls[lscnt] = c;
          if (++lscnt >= DGMAPL_LISTSIZE) {
#pragma omp critical
            { /* dump */
              for (j = 0; j < lscnt; j++)
                dgmapl_seti2(arr[ls[j]], val);
            }
            lscnt = 0; /* clear the buffer */
          }
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
  { /* empty the list */
    for (j = 0; j < lscnt; j++)
      dgmapl_seti2(arr[ls[j]], val);
  }
  return cnt;
}



#if DGMAPL_EXISTS
/* for n = 9, nr <= 9!/18 = 20160, |fb| <= 7! = 5040
 * so they can be contained in a 16-bit integer */
typedef struct {
  int n;
  int k;
  dgmapl_fb_t (*fbnr)[2];
} dgmapl_t;



INLINE dgmapl_t *dgmapl_open(int n, int k)
{
  dgmapl_t *mapl;
  uint64_t size;
  static const int kdef[] = {0,
    0, 1, 2, 3, 4, 4, 5, 6, 8, 9};
  clock_t t0;

  die_if (n > DGMAPL_NMAX, "n %d is too large\n", n);
  xnew(mapl, 1);
  mapl->n = n;
  if (k <= 0) /* choose a default k */
    k = kdef[n];
  mapl->k = k;
  size = (uint64_t) 1u << (n * (n - 1) / 2 - k);
  die_if (n*(n-1)/2-k > (int) sizeof(code_t) * 8,
      "n %d, k %d cannot be contained in %d bits, increase CODEBITS\n",
      n, k, (int) sizeof(code_t));
  xnew(mapl->fbnr, size);
  fprintf(stderr, "%4d: dgmapl allocated %gGB memory\n", inode,
      1. * size * sizeof(mapl->fbnr[0]) / (1024.*1024.*1024.) );
  t0 = clock();
  /* since DGMAPL_BAD is 0x8080, 0x808080, or 0x80808080, we can use this trick */
  memset(mapl->fbnr, 0x80, size * sizeof(mapl->fbnr[0]));
  fprintf(stderr, "%4d: dgmapl memory set, %gs\n",
      inode, 1.*(clock() - t0)/CLOCKS_PER_SEC);
  return mapl;
}



INLINE void dgmapl_close(dgmapl_t *mapl)
{
  free(mapl->fbnr);
  free(mapl);
}



#define dgmapl_fbnr_lookup(g, k, nr) dgmapl_fbnr_lookup0(g, k, nr, 0, NULL, NULL)

/* compute fb and nr by the larger lookup table
 * set k = 0 to use the default search length */
INLINE double dgmapl_fbnr_lookup0(const dg_t *g, int k, double *nr,
    int nocsep, int *ned, int *degs)
{
#ifdef DGMAPL_DEBUG
  static uint64_t mis = 0, tot = 0, cnt = 0, hit = 0;
#endif
  static dgmapl_t *mapl[DGMAPL_NMAX + 1];
  static int st[DGMAPL_NMAX + 1];
#pragma omp threadprivate(st)

  code_t c;
  int hasnew;
  DG_DEFN_(g);
  dgmapl_int_t ifbnr[2];
  dgmapl_t *tb;
  double fb;

  die_if (n <= 0 || n > DGMAPL_NMAX, "bad n %d\n", n);

  /* allocate memory */
  if (mapl[DG_N_] == NULL) {
#ifdef _OPENMP
    /* We test mapl[n] twice here. If mapl[n] != NULL, then the map
     * memory has been allocated, and we may safely skip the critical
     * block.  If mapl[n] == NULL, the if clause inside the critical
     * block ensures that the memory is not allocated twice.  That is,
     * after the first thread leaves the critical block, the second
     * thread can see `mapl[n]' is no longer NULL, and skip the block */
#pragma omp critical
    {
      if (mapl[DG_N_] == NULL) {
#endif /* _OPENMP */
        clock_t t0 = clock();
        mapl[DG_N_] = dgmapl_open(DG_N_, k);
        fprintf(stderr, "%4d: initialized fb/nr-mapl for %d, k %d, time %gs\n",
            inode, DG_N_, mapl[DG_N_]->k, 1.*(clock() - t0)/CLOCKS_PER_SEC);
#ifdef _OPENMP
      }
    } /* omp critical */
#endif /* _OPENMP */
  }

  tb = mapl[DG_N_];
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

  dgmapl_geti2(ifbnr, tb->fbnr[c]);

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
      /* the ring content is nonnegative */
      ifbnr[1] = (dgmapl_int_t) (*nr + .5);
      hasnew = 1;
    }
  }

  if (hasnew) {
#ifdef DGMAPL_SIMPLEUPDATE
#pragma omp critical
    {
      dgmapl_seti2(tb->fbnr[c], ifbnr);
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

#endif /* DGMAPL_EXIST */


#endif

