#ifndef DGHASH0_H__
#define DGHASH0_H__

/* randomized hash code */

#ifndef N
#error "must define N to use the hash table"
#endif

#ifndef CODEBITS
#if N >= 9
#define CODEBITS 64
#endif
#endif

#include "dg.h"
#include "dgaut.h"
#include "dgring.h"




/* the number of words to save the code of a graph of vertices
 * the number of bits needed is n * (n - 1) / 2
 * divide this by sizeof(code_t) is DGHASH_CWORDS */
#define DGHASH_CBITS  (N*(N-1)/2)
#define DGHASH_CWORDS ((DGHASH_CBITS + CODEBITS - 1)/CODEBITS)


/* linear congruential generator */
#if CODEBITS == 32
#define DGHASH_LCG(c) ((c) * 314159265u + 1u)
#else
#define DGHASH_LCG(c) ((c) * 6364136223846793005ull + 1442695040888963407ull)
#endif


/* we will use a random number generator to scramble to code */
#if DGHASH_CWORDS == 1
  #define DGHASH_GETHASHID(id, c, hbits) \
    id = (DGHASH_LCG(c[0]) >> (CODEBITS - hbits))
#else
  #define DGHASH_GETHASHID(id, c, hbits) { int k_; \
    for (id = DGHASH_LCG(c[0]), k_ = 1; k_ < DGHASH_CWORDS; k_++) \
      id = DGHASH_LCG(id + c[k_]); \
    id >>= CODEBITS - hbits; }
#endif /* DGHASH_CWORDS == 1 */



/* copy the code */
#if DGHASH_CWORDS == 1
  #define DGLS_CCPY(dest, src) dest[0] = src[0]
#elif DGHASH_CWORDS == 2
  #define DGLS_CCPY(dest, src) { \
    dest[0] = src[0]; dest[1] = src[1]; }
#else
  #define DGLS_CCPY(dest, src) { int cpk_; \
    for (cpk_ = 0; cpk_ < DGHASH_CWORDS; cpk_++) \
      dest[cpk_] = src[cpk_]; }
#endif



/* compare two scodes */
INLINE int dgls_ccmp(const code_t *a, const code_t *b)
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
  for (k = DGHASH_CWORDS - 1; k >= 0; k--) {
    if (a[k] > b[k]) return 1;
    else if (a[k] < b[k]) return -1;
  }
#endif
  return 0;
}



#if N <= 13
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
  dgls_fb_t fb;
  dgls_fb_t nr;
} dglsent_t;

typedef struct {
  int cnt;
  int cap;
  dglsent_t *arr;
} dgls_t;

typedef struct {
  int n;
  int bits; /* number of hash bits */
  size_t blksz;
  size_t lsn; /* number of lists */
  dgls_t *ls;
  size_t mem; /* total memory usage in bytes */
  size_t memmax; /* maximal memory usage in bytes */
} dghash_t;



INLINE dghash_t *dghash_open(int n, int bits, size_t blksz, size_t memmax)
{
  /* default parameters */
  static struct { int bits, blksz; } defp[] = {{1, 1},
    {1, 1}, {1, 1}, {1, 1}, {3, 2}, {5, 2},
    {6, 2}, {6, 4}, {10, 4}, {14, 4}, {20, 4},
  };
  dghash_t *h;
  size_t i;

  /* for n > 10, we use 20 bit, which gives a relatively small hash table
   * this allows a more efficient use of limited memory, because all
   * hash table entries will be used */
  if (bits <= 0) bits = (n <= 10) ? defp[n].bits : 20;
  if (blksz == 0) blksz = (n <= 10) ? defp[n].blksz : 4;
  if (memmax == 0) memmax = 0x40000000; /* 1GB */
  fprintf(stderr, "%4d: bits %d, blksz %u, mem %gM, cbits %d, cwords %d\n",
      inode, bits, (unsigned) blksz, memmax/(1024.*1024), DGHASH_CBITS, DGHASH_CWORDS);

#ifdef N
  die_if (n != N, "n %d != N %d\n", n, N);
#endif
  xnew(h, 1);
  h->n = n;
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
INLINE code_t dghash_getid(const dg_t *g, code_t *c, int hbits)
{
  code_t id;
  static dg_t *ng;
#pragma omp threadprivate(ng);

  if (ng == NULL) { /* allocate ng */
    DG_DEFN_(g);
#ifdef N
    ng = dg_open(N);
#else
    ng = dg_open(DG_NMAX);
#endif
  }
  dg_canlabel(ng, g);
  dg_encode(ng, c);
  DGHASH_GETHASHID(id, c, hbits);
  return id;
}


/* find the code_t *c in the list
 * if it fails, *pos returns the point for insertion */
INLINE int dgls_find(dgls_t *ls, code_t *c, int *pos)
{
  int imin, imax, imid, cmp;

  /* start the binary search */
  imin = 0;
  imax = ls->cnt - 1;
  while (imax >= imin) {
    imid = (imin + imax) / 2;
    cmp = dgls_ccmp(c, ls->arr[imid].c);
    if (cmp == 0) return imid;
    if (cmp > 0) imin = imid + 1;
    else imax = imid - 1;
  }
  *pos = imin;
  return -1;
}



/* insert the new entry just before position ipos */
INLINE int dgls_add(dgls_t *ls, int ipos, const code_t *c,
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
  DGLS_CCPY(ent->c, c);
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
  static dghash_t *hash[DG_NMAX + 1];
  code_t c[DGHASH_CWORDS];
#pragma omp threadprivate(hash, c)
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

  ls = h->ls + dghash_getid(g, c, h->bits);
  pos = dgls_find(ls, c, &ipos);
  if (pos >= 0) { /* entry exists */
    *nr = ls->arr[pos].nr;
    return ls->arr[pos].fb;
  }
  fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
  *nr = dg_nring_mixed0(g, ned, degs);
#pragma omp critical
  {
    dgls_add(ls, ipos, c, fb, *nr, h->blksz, &h->mem, h->memmax);
  }
  return fb;
}



/* compute and print statistics of the usage of the hash table */
INLINE void dghash_stat(const dghash_t *h, FILE *fp)
{
  size_t i, nz = 0;
  int x, hmax = 0, hmin = 10000;
  double cnt = 0, sm = 0, sm2 = 0, smb = 0, sm2b = 0;

  for (i = 0; i < h->lsn; i++) {
    x = h->ls[i].cnt;
    sm += x;
    sm2 += x * x;
    if (x > 0) {
      cnt += x;
      nz += 1;
      smb += x;
      sm2b += x * x;
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
  fprintf(fp, "cnt %.0f, ls: %" PRIu64 "/%" PRIu64 " (%5.2f%%), "
      "av. %.3f(%.3f), nzav. %.2f(%.2f), %d-%d, mem. %.2fM\n",
      cnt, (uint64_t) nz, (uint64_t) h->lsn, 100.*nz/h->lsn,
      sm, sqrt(sm2), smb, sqrt(sm2b), hmin, hmax, h->mem/(1024.*1024));
}



#endif /* DGHASH_H__ */

