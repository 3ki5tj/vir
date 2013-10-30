#ifndef DGHASH_H__
#define DGHASH_H__

/* full hash table */

#ifndef DGHASH_EXISTS
#define DGHASH_EXISTS 1
#endif


#ifdef DGHASH_EXISTS


/* although set code bits to 64 would possible make hash function faster
 * it might make Wheatley's method slow, esp on a 32-bit machine
 * the disadvantage outweights the advantage */
/*
#ifndef CODEBITS
#if defined(N) && (N >= 9)
#define CODEBITS 64
#endif
#endif
*/

#include "dg.h"
#include "dgaut.h"
#include "dgring.h"




/* the number of words to save the code of a graph of vertices
 * the number of bits needed is n * (n - 1) / 2
 * divide this by sizeof(code_t) is DGHASH_CWORDS */
#define DGHASH_CBITS  (DG_NMAX *(DG_NMAX - 1)/2)
#define DGHASH_CWORDS ((DGHASH_CBITS + CODEBITS - 1)/CODEBITS)


/* linear congruential generator */
#if CODEBITS == 32
#define DGHASH_LCG(c) ((c) * 314159265u + 1u)
#else
#define DGHASH_LCG(c) ((c) * 6364136223846793005ull + 1442695040888963407ull)
#endif


/* we will use a random number generator to scramble to code */
#if DGHASH_CWORDS == 1
  #define DGHASH_GETID(id, c, cnt, hbits) \
    id = (DGHASH_LCG(c[0]) >> (CODEBITS - hbits))
#elif DGHASH_CWORDS == 2
  #define DGHASH_GETID(id, c, cnt, hbits) \
    id = DGHASH_LCG(DGHASH_LCG(c[0]) + c[1]) >> (CODEBITS - hbits)
#else
  #define DGHASH_GETID(id, c, cnt, hbits) { int k_; \
    for (id = DGHASH_LCG(c[0]), k_ = 1; k_ < cnt; k_++) \
      id = DGHASH_LCG(id + c[k_]); \
    id >>= CODEBITS - hbits; }
#endif /* DGHASH_CWORDS == 1 */

/* map DGHASH_CWORDS_ to the macro if possible, or a variable otherwise */
#ifdef N
#define DGHASH_CWORDS_ DGHASH_CWORDS
#else
#define DGHASH_CWORDS_ cwords
#endif



/* copy the code */
#if DGHASH_CWORDS == 1
  #define DGLS_CCPY(dest, src, cnt) dest[0] = src[0]
#elif DGHASH_CWORDS == 2
  #define DGLS_CCPY(dest, src, cnt) { \
    dest[0] = src[0]; dest[1] = src[1]; }
#else
  #define DGLS_CCPY(dest, src, cnt) { int cpk_; \
    for (cpk_ = 0; cpk_ < cnt; cpk_++) \
      dest[cpk_] = src[cpk_]; }
#endif



/* compare two scodes */
INLINE int dgls_ccmp(const code_t *a, const code_t *b, int cwords)
{
#if DGHASH_CWORDS == 1
  if (a[0] > b[0]) return 1;
  else if (a[0] < b[0]) return -1;
#elif DGHASH_CWORDS == 2
  if (a[1] > b[1]) return 1;
  if (a[1] < b[1]) return -1;
  if (a[0] > b[0]) return 1;
  if (a[0] < b[0]) return -1;
#else
  int k;
  for (k = DGHASH_CWORDS_ - 1; k >= 0; k--) {
    if (a[k] > b[k]) return 1;
    else if (a[k] < b[k]) return -1;
  }
#endif
  return 0;
}



#if defined(N) && (N <= 13)
  /* this saves some storage space */
  typedef int32_t dgls_fb_t;
  /* convert double to integer and save */
  #define DGLS_SAVEFB(a, b) \
    (a) = (dgls_fb_t) (((b) < 0 ? ((b) - .5) : ((b) + .5)))
#else
  typedef double dgls_fb_t;
  #define DGLS_SAVEFB(a, b) (a) = (b)
#endif


typedef struct {
  code_t c[DGHASH_CWORDS];
  dgls_fb_t fb, nr;
} dglsent_t;

typedef struct {
  int cnt;
  int cap;
  dglsent_t *arr;
} dgls_t;

typedef struct {
  int n;
  int bits; /* number of hash bits */
  int cwords; /* number of code_t to save the connectivity matrix */
  size_t blksz;
  size_t lsn; /* number of lists */
  dgls_t *ls;
  size_t mem; /* total memory usage in bytes */
  size_t memmax; /* maximal memory usage in bytes */
  int dostat;
  double tot, hits; /* statistics */
} dghash_t;



INLINE dghash_t *dghash_open(int n, int bits, size_t blksz, size_t memmax)
{
  /* default parameters */
  static struct { int bits, blksz; } defp[] = {{1, 1},
    {1, 1}, {1, 1}, {1, 1}, {3, 2}, {5, 2},
    {6, 2}, {6, 4}, {10, 4}, {12, 4}, {18, 4},
  };
  dghash_t *h;
  size_t i;

  /* for n > 10, we use 20 bit, which gives a relatively small hash table
   * this allows a more efficient use of limited memory, because all
   * hash table entries will be used */
  if (bits <= 0) bits = (n <= 10) ? defp[n].bits : 20;
  if (blksz == 0) blksz = (n <= 10) ? defp[n].blksz : 4;
  if (memmax == 0) memmax = 0x40000000; /* 1GB */

  xnew(h, 1);
  h->n = n;
  h->cwords = (n * (n - 1) / 2 + CODEBITS - 1) / CODEBITS;
  h->bits = bits;
  h->blksz = blksz;
  h->memmax = memmax;
  h->lsn = (size_t) 1u << bits;
  xnew(h->ls, h->lsn);
  for (i = 0; i < h->lsn; i++) {
    h->ls[i].cnt = h->ls[i].cap = 0;
    h->ls[i].arr = NULL;
  }
  h->mem = h->lsn * sizeof(h->ls[0]);
  if (inode == MASTER) {
    fprintf(stderr, "dghash bits %d, blksz %u, mem %gM, cbits %d, cwords %d, unit %dB\n",
        bits, (unsigned) blksz, memmax/(1024.*1024),
        n*(n - 1)/2, h->cwords, (int) sizeof(dglsent_t));
    if (h->cwords * 2 < DGHASH_CWORDS)
      fprintf(stderr, "dghash can save memory with predefined N (%d >> %d) or DGHASH_CWORDS (%d >> %d)\n",
          DG_NMAX, n, DGHASH_CWORDS, h->cwords);
  }
  h->dostat = 0; /* turn off stat by default, the user can turn it on manually */
  h->tot = h->hits = 1e-30;
  return h;
}



INLINE void dghash_close(dghash_t *h)
{
  size_t i;

  for (i = 0; i < h->lsn; i++)
    if (h->ls[i].arr != NULL)
      free(h->ls[i].arr);
  free(h->ls);
  free(h);
}



/* return the hash list id */
INLINE code_t dghash_getid(const dg_t *g, code_t *c, int cwords, int hbits)
{
  code_t id;
  static dg_t *ng;
#pragma omp threadprivate(ng)

  if (ng == NULL) { /* allocate ng */
    ng = dg_open(DG_NMAX);
  }
  dg_canlabel(ng, g);
  dg_encode(ng, c);
  DGHASH_GETID(id, c, DGHASH_CWORDS_, hbits);
  return id;
}


/* find the code_t *c in the list
 * if it fails, *pos returns the point for insertion */
INLINE int dgls_find(dgls_t *ls, code_t *c, int cwords, int *pos)
{
  int imin, imax, imid, cmp;

  /* start the binary search */
  imin = 0;
  imax = ls->cnt - 1;
  while (imax >= imin) {
    imid = (imin + imax) / 2;
    cmp = dgls_ccmp(c, ls->arr[imid].c, DGHASH_CWORDS_);
    if (cmp == 0) return imid;
    if (cmp > 0) imin = imid + 1;
    else imax = imid - 1;
  }
  *pos = imin;
  return -1;
}



/* insert the new entry just before position ipos */
INLINE int dgls_add(dgls_t *ls, int ipos, const code_t *c, int cwords,
    double fb, double nr, size_t blksz, size_t *mem, size_t memmax)
{
  dglsent_t *ent;
  if (ls->cap == 0) {
    if (*mem > memmax) return -1;
    *mem += blksz * sizeof(dglsent_t);
    ls->cap = blksz;
    xnew(ls->arr, ls->cap);
  } else if (ls->cnt >= ls->cap) { /* allocate memory */
    if (*mem > memmax) return -1;
    *mem += blksz * sizeof(dglsent_t);
    ls->cap += blksz;
    xrenew(ls->arr, ls->cap);
  }
  ent = ls->arr + ipos;
  if (ipos < ls->cnt)
    memmove(ent + 1, ent, (ls->cnt - ipos) * sizeof(dglsent_t));
  DGLS_CCPY(ent->c, c, DGHASH_CWORDS_);
  DGLS_SAVEFB(ent->fb, fb);
  DGLS_SAVEFB(ent->nr, nr);
  ls->cnt++;
  return 0;
}


#define dghash_fbnr_lookup(g, nr) \
  dghash_fbnr_lookup0(NULL, g, nr, 0, NULL, NULL)

INLINE double dghash_fbnr_lookup0(dghash_t *h, const dg_t *g,
    double *nr, int nocsep, int *ned, int *degs)
{
  static code_t c[DGHASH_CWORDS];
#pragma omp threadprivate(c)
  static dghash_t *hash[DG_NMAX + 1]; /* default hash, shared */
  DG_DEFN_(g);
  int pos, ipos;
  dgls_t *ls;
  double fb;

  if (h == NULL) { /* initialize the stock hash table */
    die_if (DG_N_ > DG_NMAX, "n %d is too large %d\n", DG_N_, DG_NMAX);
    if (hash[DG_N_] == NULL) {
#ifdef _OPENMP
#pragma omp critical
      {
        if (hash[DG_N_] == NULL)
#endif /* _OPENMP */
          hash[DG_N_] = dghash_open(DG_N_, 0, 0, 0);
#ifdef _OPENMP
      } /* omp critical */
#endif /* _OPENMP */
    }
    h = hash[DG_N_];
  }

  ls = h->ls + dghash_getid(g, c, h->cwords, h->bits);
  pos = dgls_find(ls, c, h->cwords, &ipos);
  if (h->dostat) { /* accumulate statistics */
#pragma omp critical
    {
      h->tot += 1;
      h->hits += (pos >= 0);
    }
  }
  if (pos >= 0) { /* entry exists */
    *nr = ls->arr[pos].nr;
    return ls->arr[pos].fb;
  }
  fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
  *nr = dg_nring_mixed0(g, ned, degs);
#pragma omp critical
  {
    dgls_add(ls, ipos, c, h->cwords,
        fb, *nr, h->blksz, &h->mem, h->memmax);
  }
  return fb;
}



/* compute and print statistics of the usage of the hash table */
INLINE void dghash_printstat(dghash_t *h, FILE *fp)
{
  size_t i, nz = 0;
  int x, hmax = 0, hmin = 10000;
  double cnt = 0, cap = 0, sm = 0, sm2 = 0, smb = 0, sm2b = 0;

  for (i = 0; i < h->lsn; i++) {
    x = h->ls[i].cnt;
    cap += h->ls[i].cap;
    sm += x;
    sm2 += 1. * x * x;
    if (x > 0) {
      cnt += x;
      nz += 1;
      smb += x;
      sm2b += 1. * x * x;
      if (x > hmax) hmax = x;
      if (x < hmin) hmin = x;
    }
  }
  sm /= h->lsn;
  sm2 = (sm2 / h->lsn) - sm * sm;
  if (nz) {
    smb /= nz;
    sm2b = (sm2b / nz) - smb * smb;
  }
  fprintf(fp, "%4d: dghash: ", inode);
  if (h->dostat) {
    fprintf(fp, "hits %5.2f%%, ", 100.*h->hits/h->tot);
#pragma omp critical
    { /* reset hits */
      h->hits = h->tot = 1e-30;
    }
  }
  fprintf(fp, "cnt %.0f, used %" PRIu64 "/%" PRIu64 " (%5.2f%%), "
      "av. %.3f(%.3f), nzav. %.2f(%.2f), %d-%d, mem. %.2fM (pack %5.2f%%)\n",
      cnt, (uint64_t) nz, (uint64_t) h->lsn, 100.*nz/h->lsn,
      sm, sqrt(sm2), smb, sqrt(sm2b), hmin, hmax,
      h->mem/(1024.*1024), 100*cnt/cap);
}

#endif /* defined(DGHASH_EXISTS) */

#endif /* DGHASH_H__ */

