#ifndef DGHASH_H__
#define DGHASH_H__

/* full hash table */


/* hash table is pretty much universal,
 * so we enable it whenever this file is included */
#ifndef DGHASH_EXISTS
#define DGHASH_EXISTS 1
#endif

#ifndef DGHASH_NMAX
#define DGHASH_NMAX DG_NMAX
#endif


#ifdef DGHASH_EXISTS


/* although set code bits to 64 would possible make hash function faster
 * it might make Wheatley's method slow on a 32-bit machine
 * the disadvantage outweights the advantage */
/*
#ifndef DG_WORDBITS
#if defined(N) && (N >= 9)
#define DG_WORDBITS 64
#endif
#endif
*/

#define ZCOM_PICK
#define ZCOM_BLKMEM /* block memory allocator */
#include "zcom.h"
#include "dg.h"
#include "dgaut.h"
#include "dgring.h"




/* the number of words to save the code of a graph of vertices
 * the number of bits needed is n * (n - 1) / 2
 * divide this by sizeof(dgword_t) is DGHASH_CWORDS */
#define DGHASH_CBITS  (DGHASH_NMAX *(DGHASH_NMAX - 1)/2)
#define DGHASH_CWORDS ((DGHASH_CBITS + DG_WORDBITS - 1)/DG_WORDBITS)


/* linear congruential generator */
#if DG_WORDBITS == 32
#define DGHASH_LCG(c) ((c) * 314159265u + 1u)
#else
#define DGHASH_LCG(c) ((c) * CU64(6364136223846793005) + CU64(1442695040888963407))
#endif


/* we will use a random number generator to scramble to code */
#if DGHASH_CWORDS == 1
  #define DGHASH_GETID(id, c, cnt, hbits) \
    id = (DGHASH_LCG(c[0]) >> (DG_WORDBITS - hbits))
#elif DGHASH_CWORDS == 2
  #define DGHASH_GETID(id, c, cnt, hbits) \
    id = DGHASH_LCG(DGHASH_LCG(c[0]) + c[1]) >> (DG_WORDBITS - hbits)
#else
  #define DGHASH_GETID(id, c, cnt, hbits) { int k_; \
    for (id = DGHASH_LCG(c[0]), k_ = 1; k_ < cnt; k_++) \
      id = DGHASH_LCG(id + c[k_]); \
    id >>= DG_WORDBITS - hbits; }
#endif /* DGHASH_CWORDS == 1 */

/* map DGHASH_CWORDS_ to the macro if possible, or a variable otherwise */
#ifdef N
#define DGHASH_CWORDS_(cwords) DGHASH_CWORDS
#else
#define DGHASH_CWORDS_(cwords) cwords
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



/* the maximal fb is (n - 2)!, and the maximal nr is n!/(2 n) > (n - 2)!
 * a 32-bit signed integer can hold a number up to 2.1e9
 * a 32-bit unsigned integer can hold a number up to 4.2e9
 * the cutoff N is 14 because 13!/(2*13) = 2.4e8, 14!/(2*14) = 3.1e9
 * using 32-bit integers can save memory */
#if defined(N) && (N <= 14)
  typedef int32_t dgls_fb_t;
  typedef uint32_t dgls_nr_t;
  /* convert double to integer and save */
  #define DGLS_SAVEFB(a, b) \
    (a) = (dgls_fb_t) (((b) < 0 ? ((b) - .5) : ((b) + .5)))
  #define DGLS_SAVENR(a, b) \
    (a) = (dgls_nr_t) ((b) + .5)
#else
  typedef double dgls_fb_t;
  typedef double dgls_nr_t;
  #define DGLS_SAVEFB(a, b) (a) = (b)
  #define DGLS_SAVENR(a, b) (a) = (b)
#endif



typedef struct {
  dgword_t c[DGHASH_CWORDS];
  dgls_fb_t fb;
#ifndef DG_NORING
  dgls_nr_t nr;
#endif
} dgitem_t;

typedef struct tagdgls_t {
  int cnt;
  dgitem_t *arr;
  struct tagdgls_t *next;
} dgls_t;

typedef struct {
  int n; /* number of vertices in the graph */
  int bits; /* number of hash bits */
  int cwords; /* number of dgword_t to save the connectivity matrix */
  size_t blksz; /* number of items to allocate for each list */
  size_t lsn; /* number of lists */
  dgls_t *ls; /* each ls[key], where key = 0..lsn-1, is a list of items
                 under the same hash key */
  size_t mem; /* total memory usage in bytes */
  size_t memmax; /* maximal memory usage in bytes */
  int isoenum;  /* enumerate other isomorphic graphs after a new graph is found */
  int isomax; /* maximal number graphs used in the above process */
  int level; /* automorphism level, see dgaut.h/dg_getrep() for details
                0 means no tranformation to the graph
                1 means to find a version compatible with the degree sequence
                4 means to find the canonical label */
  int dostat; /* collect the following data */
  double tot, hits; /* statistics */
  /* memory allocators, since memory allocation are done
   * in omp critical sections, we don't need to make them thread private */
  blkmem_t *blkmem_dgls, *blkmem_dgitem;
} dghash_t;




/* open an hash table for the nth virial coefficient
 * `bits'     is the number of bits of the hash table, such that it has
 *            2^bits different keys, hence buckets
 *            each bucket poccess a linked list
 *            0 or negative values: for the default value
 * `blksz'    is the number of blocks in a link of each hash table list
 *            0: for the default value
 * `blkmem'   is the number of memory blocks used in the underlying
 *            pooled memory allocator for the linked lists
 *            0: for the default value
 * `memmax'   is the maximal permissble memory in bytes. Even if it is infinity
 *            some latter allocation through the blkmem_new() function may fail
 *            but the program should not abort
 *            However, using too much memory may slow down the program
 *            0: for the default value
 * `initls'   means to try to initially allocate spaces for all keys
 *            0: off, 1: on, -1: default
 * `level'    means the type of isomorphic transformation to the graph
 *            before using the hash table
 *            possible values 0-4; a value > 9 or a negative value means
 *            to use the canonical label, the value 1 means to use a
 *            permutation compatible with the degree sequence
 * `isoenum'  0 or 1, means to find others ways of permutating vertices
 *            it is unnecessary if the canonical label is used (level == 4)
 *            but useful if the degree sequence is used (level == 1)
 * `isomax'   is the maximal number of isomorphic graphs in the above search
 *            for each newly discovered graph */
INLINE dghash_t *dghash_open(int n, int bits,
    size_t blksz, size_t blkmem, size_t memmax, int initls,
    int level, int isoenum, int isomax)
{
  /* generally blksz should be large enough to avoid too many linked lists */
  /* default parameters for level 1 */
  static struct { int bits, blksz; } defp1[] = {
    {26, 4}, /* default setting for n > 10 */
    { 1, 1}, { 1, 1}, { 1, 1}, { 3, 2}, { 6, 2},
    {10, 4}, {14, 4}, {20, 4}, {24, 4}, {26, 4},
  };
  /* default parameters for level 4 */
  static struct { int bits, blksz; } defp4[] = {
    {24, 4}, /* default setting for n > 10 */
    { 1, 1}, { 1, 1}, { 1, 1}, { 3, 2}, { 6, 2},
    { 8, 2}, {10, 4}, {14, 4}, {22, 4}, {24, 4},
  };
  dghash_t *h;
  size_t i;
  /* approximate hash table size (with low level of automorphism)
   *
   * level 1 uses most memory
   *
   * D        N   level  size (Mb)
   * 2        9     1      160M
   * 3        9     1      500M
   * >=4      9     1      900M
   *
   * 2       10     1    ~2000M
   * 3       10     1    ~8000M
   * 4       10     1   >13000M
   *
   * level 4 uses least memory
   *
   * 3        9     4       29M
   * 3       10     4      180M
   * */

  if (level >= 4 || level < 0) {
    if (bits    <=  0)  bits    = (n <= 10) ? defp4[n].bits    : defp4[0].bits;
    if (blksz   ==  0)  blksz   = (n <= 10) ? defp4[n].blksz   : defp4[0].blksz;
    isoenum = 0;
  } else if (level == 1) {
    if (bits    <=  0)  bits    = (n <= 10) ? defp1[n].bits    : defp1[0].bits;
    if (blksz   ==  0)  blksz   = (n <= 10) ? defp1[n].blksz   : defp1[0].blksz;
    if (isoenum <   0)  isoenum = 1;
  }
  die_if (bits > DG_WORDBITS, "dghash: too many hash bits %d > %d\n",
      bits, DG_WORDBITS);
  if (memmax  ==  0)  memmax  = 0x60000000; /* 1.5GB */
  if (isomax  <=  0)  isomax  = 1000;

  xnew(h, 1);
  h->n = n;
  h->cwords = (n * (n - 1) / 2 + DG_WORDBITS - 1) / DG_WORDBITS;
  h->bits = bits;
  h->blksz = blksz;
  h->memmax = memmax;
  h->level = level;
  h->isoenum = isoenum;
  h->isomax = isomax;

  /* blkmem is the number of objects (list or list items) we allocate in each time
   * if blkmem == 1, the system will leave many holes in the memory space, making
   * it inefficient */
  if (blkmem  ==  0)  blkmem  = 0x100; /* 65536 */
  h->blkmem_dgls = blkmem_open(sizeof(dgls_t), blkmem);
  h->blkmem_dgitem = blkmem_open(sizeof(dgitem_t), h->blksz * blkmem);

  h->lsn = (size_t) 1u << (size_t) bits;
  xnew(h->ls, h->lsn);
  for (i = 0; i < h->lsn; i++) {
    h->ls[i].cnt = 0;
    h->ls[i].next = NULL;
    h->ls[i].arr = NULL;
  }
  h->mem = sizeof(*h) + h->lsn * sizeof(h->ls[0]);

  /* Assuming hashbits are properly given, it might be advantageous
   * to allocate memory for all lists to preserve space for unused entries
   * This can save much time for latter allocations
   * This may or may not allow a better coverage in the long run */
  if (initls < 0) /* turn it on by default */
    initls = 1;
  if (h->mem + 1. * h->blksz * sizeof(dgitem_t) * h->lsn > h->memmax)
    initls = 0;
  if (initls) {
    size_t sz = h->blksz * sizeof(dgitem_t) * h->lsn;
    fprintf(stderr, "dghash: initializing all lists, %.0f items %gM/%gM\n",
        1.*h->lsn, sz / (1024.*1024), h->memmax / (1024.*1024));
    if ((h->ls[0].arr = calloc(1, sz)) != NULL) {
      for (i = 1; i < h->lsn; i++)
        h->ls[i].arr = h->ls[0].arr + i * h->blksz;
      h->mem += sz;
    } else {
      fprintf(stderr, "dghah: initial list allocation failed\n");
    }
  }

  fprintf(stderr, "dghash bits %d, blksz %u, mem %gM/%gM, "
      "cbits %d, cwords %d, unit %d bytes, level %d, isoenum %d, isomax %d\n",
      bits, (unsigned) blksz, h->mem/(1024.*1024), h->memmax/(1024.*1024),
      n*(n - 1)/2, h->cwords, (int) sizeof(dgitem_t), level, isoenum, isomax);
  if (h->cwords * 2 < DGHASH_CWORDS)
    fprintf(stderr, "dghash can save memory with predefined N (%d >> %d) or DGHASH_CWORDS (%d >> %d)\n",
        DGHASH_NMAX, n, DGHASH_CWORDS, h->cwords);

  h->dostat = 1;
#ifdef _OPENMP
  /* statistics greatly affects the performance,
   * turn off stat by default, the user can turn it on manually */
  if (nnodes > 1) h->dostat = 0;
#endif
  h->tot = h->hits = 1e-30;
  return h;
}



/* in thread case this function should not be called for safety */
INLINE void dghash_close(dghash_t *h)
{
#if 0 /* to print allocation statistics */
  blkmem_print(h->blkmem_dgitem, "dgitem");
  blkmem_print(h->blkmem_dgls, "dgls");
#endif
  blkmem_close(h->blkmem_dgitem);
  blkmem_close(h->blkmem_dgls);
  free(h->ls);
  free(h);
}



/* find the dgword_t *c in the list
 * This function do not have to be placed an omp critical block
 * This, however, requires that we add items sequentially
 * and we will use linear search instead of binary search */
INLINE dgls_t *dgls_find(dgls_t *ls, dgword_t *c, int cn, int *pos)
{
  int i, lscnt;

  /* loop over the linked list */
  while ( ls != NULL ) {
    lscnt = ls->cnt;
    for (i = 0; i < lscnt; i++) {
      //die_if (ls->arr[i].c[0] == 0, "uninitialied array element i %d\n", i);
      if ( dgcode_eq(c, ls->arr[i].c, cn) ) {
        *pos = i;
        return ls;
      }
    }
    ls = ls->next;
  }
  return NULL;
}



/* add the new entry the list `ls'
 * this function is expected to be included in a omp critical block */
INLINE int dgls_add(dgls_t *ls, const dgword_t *c,
    double fb, double nr, dghash_t *h)
{
  dgitem_t *item;
  int lscnt;

  /* allocate memory if needed */
  /* loop till the end of the linked list */
  while (ls->next != NULL)
    ls = ls->next;
  lscnt = ls->cnt;
  //die_if (lscnt > 0 && lscnt > blksz, "cnt %d, cap %d\n", lscnt, blksz);
  if (ls->arr == 0) { /* no item in the list */
    if (h->mem > h->memmax) return -1;
    h->mem += h->blksz * sizeof(ls->arr[0]);
    if ((ls->arr = blkmem_new(h->blkmem_dgitem, h->blksz)) == NULL)
      return -1;
    /* this array will not be reallocated */
    //xnew(ls->arr, h->blksz);
    //lscnt = 0; /* to be incremented below */
  } else {
    /* check if we have exhausted the capacity of this list */
    if ((size_t) lscnt == h->blksz) {
      if (h->mem > h->memmax) return -1;
      h->mem += h->blksz * sizeof(ls->arr[0]) + sizeof(*ls);
      //xnew(ls->next, 1);
      if ((ls->next = blkmem_new(h->blkmem_dgls, 1)) == NULL)
        return -1;
      //xnew(ls->next->arr, h->blksz);
      if ((ls->next->arr = blkmem_new(h->blkmem_dgitem, h->blksz)) == NULL) {
        ls->next = NULL;
        return -1;
      }
      ls = ls->next; /* shift to the next list */
      lscnt = 0; /* to be incremented below */
    }
  }

  item = ls->arr + lscnt;
  DGCODE_CPY(item->c, c, h->cwords);
  DGLS_SAVEFB(item->fb, fb);
#ifndef DG_NORING
  DGLS_SAVENR(item->nr, nr);
#endif
  ls->cnt = lscnt + 1;
#pragma omp flush
  return 0;
}



static int dghash_isocap_;
static dgword_t *dghash_isocodes_;
#pragma omp threadprivate(dghash_isocodes_, dghash_isocap_)

/* add other keys compatbile with the automorphism level into the hash list */
INLINE int dghash_enumiso(dghash_t *h, const dg_t *g,
    double fb, double nr, dg_t *ng)
{
  int i, cnt, icnt, pos;
  int cwords = h->cwords, hbits = h->bits;
  dgword_t id, *c;
  dgls_t *ls, *ls1;

  if (cwords * h->isomax > dghash_isocap_) {
    if (dghash_isocodes_) free(dghash_isocodes_);
    dghash_isocap_ = cwords * h->isomax;
    xnew(dghash_isocodes_, dghash_isocap_);
#pragma omp critical
    {
      h->mem += dghash_isocap_ * sizeof(dgword_t);
    }
  }
  cnt = dg_repisocodels(dghash_isocodes_, cwords, h->isomax, g, h->level, ng);
#pragma omp critical
  {
    for (icnt = 0, i = 0; i < cnt; i++) {
      c = dghash_isocodes_ + i * cwords;
      DGHASH_GETID(id, c, cwords, hbits);
      ls = h->ls + id;
      if ((ls1 = dgls_find(ls, c, cwords, &pos)) != NULL) /* item exists */
        continue;
      dgls_add(ls, c, fb, nr, h);
      icnt++;
    }
  }
  return icnt; /* number of items added to the list */
}



#define dghash_fbnr_lookup(g, nr) \
  dghash_fbnr_lookup0(NULL, g, nr, 0, NULL, NULL)

static dghash_t *dghash_[DGHASH_NMAX + 1]; /* default hash, shared */

INLINE double dghash_fbnr_lookup0(dghash_t *h, const dg_t *g,
    double *nr, int nocsep, int *ned, int *degs)
{
  static dgword_t c[DGHASH_CWORDS];
  static dgword_t ng_c[DG_NMAX];
  static dg_t ng[1] = {{DGHASH_NMAX, NULL}}; /* a stock graph */
#pragma omp threadprivate(c, ng_c, ng)

  DG_DEFN_(g)
  int pos;
  dgword_t hashid;
  dgls_t *ls, *ls1;
  double fb = 0;

  if (h == NULL) { /* initialize the stock hash table */
    die_if (DG_N_ > DGHASH_NMAX, "n %d is too large %d\n", DG_N_, DGHASH_NMAX);
    if (dghash_[DG_N_] == NULL) {
#ifdef _OPENMP
#pragma omp critical
      {
        if (dghash_[DG_N_] == NULL)
#endif /* _OPENMP */
          /* open a hash table with default settings */
          dghash_[DG_N_] = dghash_open(DG_N_, 0, 0, 0, 0, -1, -1, -1, 0);
#ifdef _OPENMP
      } /* omp critical */
#endif /* _OPENMP */
    }
    h = dghash_[DG_N_];
  }

  /* find the representative graph of according to the automorphism level */
  if (ng->c == NULL) ng->c = ng_c; /* initialize the stock graph */
  dg_repiso(ng, g, h->level);
  dg_encode(ng, c);
  DGHASH_GETID(hashid, c, DGHASH_CWORDS_(h->cwords), h->bits);
  ls = h->ls + hashid;
//#pragma omp flush /* flush to get the recent view */
  ls1 = dgls_find(ls, c, h->cwords, &pos);
  if (ls1 != NULL) { /* entry exists */
#ifndef DG_NORING
    *nr = ls1->arr[pos].nr;
#else
    *nr = 0;
#endif
    fb  = ls1->arr[pos].fb;
  }
  if (h->dostat) { /* accumulate data for statistics */
#pragma omp critical
    {
      h->tot += 1;
      h->hits += (ls1 != NULL);
    }
  }
  if (ls1 != NULL) return fb;
  fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
#ifndef DG_NORING
  *nr = dg_nring_mixed0(g, ned, degs);
#else
  *nr = 0;
#endif
  if (h->isoenum) {
    /* omp critical region will be limited in dghash_enum() */
    dghash_enumiso(h, g, fb, *nr, ng);
  } else {
#pragma omp critical
    {
      dgls_add(ls, c, fb, *nr, h);
    }
  }
  return fb;
}



/* compute and print statistics of the usage of the hash table */
INLINE void dghash_printstat(dghash_t *h, FILE *fp)
{
  size_t i, nz = 0;
  int x, xm, hmax = 0, hmin = 10000, ip;
  double cnt = 0, cap = 1e-30, sm = 0, sm2 = 0, smb = 0, sm2b = 0;
  char sbuf[256], *p = sbuf;

  for (i = 0; i < h->lsn; i++) {
    dgls_t *ls;

    x = h->ls[i].cnt;
    xm = h->blksz;
    /* loop over linked lists */
    for (ls = h->ls[i].next; ls; ls = ls->next) {
      x += ls->cnt;
      xm += h->blksz;
    }
    cap += xm;
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
  ip = sprintf(p, "%4d: dghash: ", inode);
  p += ip;
  if (h->dostat) {
#pragma omp critical
    {
      ip = sprintf(p, "hits %5.2f%%, ", 100.*h->hits/h->tot);
      p += ip;
      /* reset hits */
      h->hits = h->tot = 1e-30;
    }
  }
  if (hmin > hmax) hmin = hmax;
  sprintf(p, "cnt %.0f, used %" PRIu64 "/%" PRIu64 " (%5.2f%%), "
      "av. %.3f(%.3f), nzav. %.2f(%.2f), %d-%d, mem. %.2fM (pack %5.2f%%)\n",
      cnt, (uint64_t) nz, (uint64_t) h->lsn, 100.*nz/h->lsn,
      sm, sqrt(sm2), smb, sqrt(sm2b), hmin, hmax,
      h->mem/(1024.*1024), 100*cnt/cap);
  fputs(sbuf, fp);
}



#endif /* defined(DGHASH_EXISTS) */



/* clean up stock objects */
INLINE void dghash_free(void)
{
#ifdef DGHASH_EXISTS
  int k;

#pragma omp critical
  {
    for (k = 0; k <= DGHASH_NMAX; k++)
      if (dghash_[k] != NULL) {
        dghash_close(dghash_[k]);
        dghash_[k] = NULL;
      }
  }

  if (dghash_isocodes_ != NULL) {
    free(dghash_isocodes_);
    dghash_isocodes_ = NULL;
    dghash_isocap_ = 0;
  }
  if (dgaut_perms_ != NULL) {
    free(dgaut_perms_);
    dgaut_perms_ = NULL;
    dgaut_nperms_ = 0;
  }
#endif /* defined(DGHASH_EXISTS) */
}



#endif /* DGHASH_H__ */

