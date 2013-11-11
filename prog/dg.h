#ifndef DG_H__
#define DG_H__
/* handling diagrams in the virial expansion
 * using bitwise operations */


#ifdef __INTEL_COMPILER
  /* 161: don't complain unknown omp pragma
   * 869: unreferenced parameter (n)
   * 981: operands evaluated in unspecified order */
  #pragma warning disable 161 869 981
#elif defined(__GNUC__) && !defined(__INTEL_COMPILER)
  /* ignore all openmp pragmas */
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
  #pragma GCC diagonstic ignored "-Wunused-parameter"
#elif defined(_MSC_VER)
  /* 4068: unknown pragma
     4146: unary minus operator applied to unsigned type */
  #pragma warning(disable : 4068 4146)
#elif defined(__BORLANDC__)
  /* 8041: negating unsigned value */
  #pragma warn -8041
#endif



#ifdef MPI
#include <mpi.h>
#define MPI_MYREAL ( (sizeof(real) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT )
MPI_Comm comm = MPI_COMM_WORLD;
#endif
#define MASTER 0
int inode = MASTER, nnodes = 1;

#pragma omp threadprivate(inode)



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RV3
#include "zcom.h"
#include <time.h>



/* DG_WORDBITS used to be called CODEBITS */
#if !defined(DG_WORDBITS) && defined(CODEBITS)
#define DG_WORDBITS CODEBITS
#endif



#ifndef DG_WORDBITS
#ifdef NMAX
#define DG_WORDBITS ((NMAX + 31)/32*32)
#else
/* Note:
 * #define DG_WORDBITS (sizeof(dgword_t) * 8)
 * doesn't work, because the compiler cannot compute sizeof() in time
 * */
#define DG_WORDBITS 32
#endif /* defined NMAX */
#endif /* defined DG_WORDBITS */



/* we only support a graph with at most DG_WORDBITS vertices */
#if defined(N) && N > DG_WORDBITS
#error only support N to DG_WORDBITS
#endif



#if DG_WORDBITS == 32
typedef uint32_t dgword_t;
typedef int32_t sdgword_t;
#elif DG_WORDBITS == 64
typedef uint64_t dgword_t;
typedef int64_t sdgword_t;
#else
#error "bad DG_WORDBITS definition"
#endif



#ifdef N
#define DG_NMAX N /* we know the size */
#else
#define DG_NMAX DG_WORDBITS
#endif



/* if the user defined a fixed N, fixed length loops can be used
 * but all diagrams must be of size n
 * otherwise we define DG_N_ as n, whatever it is */
#ifdef  N
#define DG_N_ N
#define DG_GN_(g) N
#define DG_DEFN_(g)
#define DG_MASKN_ MKBITSMASK(N)
#define DG_DEFMASKN_()
#else
#define DG_N_ n
#define DG_GN_(g) ((g)->n)
#define DG_DEFN_(g) int n = (g)->n;
#define DG_MASKN_ maskn
#define DG_DEFMASKN_() dgword_t maskn = MKBITSMASK(n);
#endif



#if DG_WORDBITS == 32

/* these masks can be generated from the helper python script mask.py */
dgword_t bitsmask_[DG_WORDBITS + 1] = {0,
        0x1,        0x3,        0x7,        0xf,
       0x1f,       0x3f,       0x7f,       0xff,
      0x1ff,      0x3ff,      0x7ff,      0xfff,
     0x1fff,     0x3fff,     0x7fff,     0xffff,
    0x1ffff,    0x3ffff,    0x7ffff,    0xfffff,
   0x1fffff,   0x3fffff,   0x7fffff,   0xffffff,
  0x1ffffff,  0x3ffffff,  0x7ffffff,  0xfffffff,
 0x1fffffff, 0x3fffffff, 0x7fffffff, 0xffffffff};

dgword_t invbitsmask_[DG_WORDBITS + 1] = {0xffffffff,
 0xfffffffe, 0xfffffffc, 0xfffffff8, 0xfffffff0,
 0xffffffe0, 0xffffffc0, 0xffffff80, 0xffffff00,
 0xfffffe00, 0xfffffc00, 0xfffff800, 0xfffff000,
 0xffffe000, 0xffffc000, 0xffff8000, 0xffff0000,
 0xfffe0000, 0xfffc0000, 0xfff80000, 0xfff00000,
 0xffe00000, 0xffc00000, 0xff800000, 0xff000000,
 0xfe000000, 0xfc000000, 0xf8000000, 0xf0000000,
 0xe0000000, 0xc0000000, 0x80000000,        0x0};

#elif DG_WORDBITS == 64

dgword_t bitsmask_[DG_WORDBITS + 1] = {0,
 CU64(0x0000000000000001), CU64(0x0000000000000003), CU64(0x0000000000000007), CU64(0x000000000000000f),
 CU64(0x000000000000001f), CU64(0x000000000000003f), CU64(0x000000000000007f), CU64(0x00000000000000ff),
 CU64(0x00000000000001ff), CU64(0x00000000000003ff), CU64(0x00000000000007ff), CU64(0x0000000000000fff),
 CU64(0x0000000000001fff), CU64(0x0000000000003fff), CU64(0x0000000000007fff), CU64(0x000000000000ffff),
 CU64(0x000000000001ffff), CU64(0x000000000003ffff), CU64(0x000000000007ffff), CU64(0x00000000000fffff),
 CU64(0x00000000001fffff), CU64(0x00000000003fffff), CU64(0x00000000007fffff), CU64(0x0000000000ffffff),
 CU64(0x0000000001ffffff), CU64(0x0000000003ffffff), CU64(0x0000000007ffffff), CU64(0x000000000fffffff),
 CU64(0x000000001fffffff), CU64(0x000000003fffffff), CU64(0x000000007fffffff), CU64(0x00000000ffffffff),
 CU64(0x00000001ffffffff), CU64(0x00000003ffffffff), CU64(0x00000007ffffffff), CU64(0x0000000fffffffff),
 CU64(0x0000001fffffffff), CU64(0x0000003fffffffff), CU64(0x0000007fffffffff), CU64(0x000000ffffffffff),
 CU64(0x000001ffffffffff), CU64(0x000003ffffffffff), CU64(0x000007ffffffffff), CU64(0x00000fffffffffff),
 CU64(0x00001fffffffffff), CU64(0x00003fffffffffff), CU64(0x00007fffffffffff), CU64(0x0000ffffffffffff),
 CU64(0x0001ffffffffffff), CU64(0x0003ffffffffffff), CU64(0x0007ffffffffffff), CU64(0x000fffffffffffff),
 CU64(0x001fffffffffffff), CU64(0x003fffffffffffff), CU64(0x007fffffffffffff), CU64(0x00ffffffffffffff),
 CU64(0x01ffffffffffffff), CU64(0x03ffffffffffffff), CU64(0x07ffffffffffffff), CU64(0x0fffffffffffffff),
 CU64(0x1fffffffffffffff), CU64(0x3fffffffffffffff), CU64(0x7fffffffffffffff), CU64(0xffffffffffffffff)};

dgword_t invbitsmask_[DG_WORDBITS + 1] = {CU64(0xffffffffffffffff),
 CU64(0xfffffffffffffffe), CU64(0xfffffffffffffffc), CU64(0xfffffffffffffff8), CU64(0xfffffffffffffff0),
 CU64(0xffffffffffffffe0), CU64(0xffffffffffffffc0), CU64(0xffffffffffffff80), CU64(0xffffffffffffff00),
 CU64(0xfffffffffffffe00), CU64(0xfffffffffffffc00), CU64(0xfffffffffffff800), CU64(0xfffffffffffff000),
 CU64(0xffffffffffffe000), CU64(0xffffffffffffc000), CU64(0xffffffffffff8000), CU64(0xffffffffffff0000),
 CU64(0xfffffffffffe0000), CU64(0xfffffffffffc0000), CU64(0xfffffffffff80000), CU64(0xfffffffffff00000),
 CU64(0xffffffffffe00000), CU64(0xffffffffffc00000), CU64(0xffffffffff800000), CU64(0xffffffffff000000),
 CU64(0xfffffffffe000000), CU64(0xfffffffffc000000), CU64(0xfffffffff8000000), CU64(0xfffffffff0000000),
 CU64(0xffffffffe0000000), CU64(0xffffffffc0000000), CU64(0xffffffff80000000), CU64(0xffffffff00000000),
 CU64(0xfffffffe00000000), CU64(0xfffffffc00000000), CU64(0xfffffff800000000), CU64(0xfffffff000000000),
 CU64(0xffffffe000000000), CU64(0xffffffc000000000), CU64(0xffffff8000000000), CU64(0xffffff0000000000),
 CU64(0xfffffe0000000000), CU64(0xfffffc0000000000), CU64(0xfffff80000000000), CU64(0xfffff00000000000),
 CU64(0xffffe00000000000), CU64(0xffffc00000000000), CU64(0xffff800000000000), CU64(0xffff000000000000),
 CU64(0xfffe000000000000), CU64(0xfffc000000000000), CU64(0xfff8000000000000), CU64(0xfff0000000000000),
 CU64(0xffe0000000000000), CU64(0xffc0000000000000), CU64(0xff80000000000000), CU64(0xff00000000000000),
 CU64(0xfe00000000000000), CU64(0xfc00000000000000), CU64(0xf800000000000000), CU64(0xf000000000000000),
 CU64(0xe000000000000000), CU64(0xc000000000000000), CU64(0x8000000000000000), CU64(0x0000000000000000)};

#endif



/* make the ith bit */
#define MKBIT(n) (((dgword_t) 1u) << (dgword_t) (n))

#if 1
/* make a mask with the lowest n bits being 1s, other bits being 0s */
#define MKBITSMASK(n) bitsmask_[n]
#define mkbitsmask(n) MKBITSMASK(n) /* for compatibility */

/* make a mask with the lowest n bits being 0s, other bits being 1s */
#define MKINVBITSMASK(n) invbitsmask_[n]

#else

/* make a mask with the lowest n bits being 1 */
INLINE dgword_t mkbitsmask(int n)
{
  if (n == DG_WORDBITS) return (dgword_t) (-1);
  else return MKBIT(n) - (dgword_t) 1;
}
#endif


/* count the number of 1 bits in x
 * http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable */
INLINE int bitcount(dgword_t x)
{
/* if the lowest two bits of x are zeros, and x has n nonzero bits,
 * then x, x+1, x+2, x+3 have respectivly n, n+1, n+1, n+2 nonzero bits */
#define B2(n) n,     n+1,     n+1,     n+2
#define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
  /* how many bits in 0..255 */
  static const unsigned char bits[256] = { B6(0), B6(1), B6(1), B6(2) };
  unsigned char *p = (unsigned char *) &x;

  /* NOTE: these shortcuts assume that the lowest bits of an integers,
   * which represent first vertices, go to the first bytes, which may be false!
   * These versions work for words of no more than N bits,
   * not necessarily of exactly N bits */
#if defined(N) && (N <= 8)
  return bits[p[0]];
#elif defined(N) && (N <= 16)
  return bits[p[0]] + bits[p[1]];
#elif defined(N) && (N <= 24)
  return bits[p[0]] + bits[p[1]] + bits[p[2]];
#elif DG_WORDBITS == 32
  return bits[p[0]] + bits[p[1]] + bits[p[2]] + bits[p[3]];
#elif DG_WORDBITS == 64
  return bits[p[0]] + bits[p[1]] + bits[p[2]] + bits[p[3]]
       + bits[p[4]] + bits[p[5]] + bits[p[6]] + bits[p[7]];
#endif
}



/* invert the lowest k bits */
INLINE dgword_t bitinvert(dgword_t x, int k)
{
  return x ^ ~(~0u << k);
}



/* reverse bits */
INLINE dgword_t bitreverse(dgword_t x)
{
  static const unsigned char bitrev[256] = {
    #define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
    #define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
    #define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
    R6(0), R6(2), R6(1), R6(3) };
  dgword_t c;
  unsigned char *p = (unsigned char *) &x;
  unsigned char *q = (unsigned char *) &c;
#if (DG_WORDBITS == 32)
  q[3] = bitrev[p[0]];
  q[2] = bitrev[p[1]];
  q[1] = bitrev[p[2]];
  q[0] = bitrev[p[3]];
#else
  q[7] = bitrev[p[0]];
  q[6] = bitrev[p[1]];
  q[5] = bitrev[p[2]];
  q[4] = bitrev[p[3]];
  q[3] = bitrev[p[4]];
  q[2] = bitrev[p[5]];
  q[1] = bitrev[p[6]];
  q[0] = bitrev[p[7]];
#endif
  return c;
}



/* ********************** BITFIRST BEGINS ****************************** */



#if (DG_WORDBITS == 32)
const int bruijn32_[32] =
  { 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};

#define BIT2ID(b) bruijn32_[((b) * 0x077CB531) >> 27]

#elif (DG_WORDBITS == 64) /* 64-bit */

const int bruijn64_[64] =
  {  0,  1,  2,  7,  3, 13,  8, 19,  4, 25, 14, 28,  9, 34, 20, 40,
     5, 17, 26, 38, 15, 46, 29, 48, 10, 31, 35, 54, 21, 50, 41, 57,
    63,  6, 12, 18, 24, 27, 33, 39, 16, 37, 45, 47, 30, 53, 49, 56,
    62, 11, 23, 32, 36, 44, 52, 55, 61, 22, 43, 51, 60, 42, 59, 58};

#define BIT2ID(b) bruijn64_[((b) * CU64(0x218A392CD3D5DBF)) >> 58]

#endif /* DG_WORDBITS == 64 */



/* macro version of bitfirstlow(); */
#define BITFIRSTLOW(id, x, b) { \
  (b) = (x) & (-(x)); \
  (id) = BIT2ID(b); }

/* find the index of the lowest 1 bit
 * on return *b is the nonzero bit
 * http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup
 * ``Using de Bruijn Sequences to Index 1 in a Computer Word''
 * by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
 * http://supertech.csail.mit.edu/papers/debruijn.pdf */
INLINE int bitfirstlow(dgword_t x, dgword_t *b)
{
  (*b) = x & (-x); /* such that only the lowest 1-bit survives */
  return BIT2ID(*b);
}



#if defined(__INTEL_COMPILER) && (DG_WORDBITS == 32)

/* directly map to the intrinsic function
 * use it only when you are sure x != 0 */
#define BITFIRSTNZ(x)  _bit_scan_forward(x)

#define BITFIRST(id, x) id = (x) ? _bit_scan_forward(x) : 0;

#define bitfirst(x) ((x) ? _bit_scan_forward(x) : 0)

#else /* defined(__INTEL_COMPILER) && (DG_WORDBITS == 32) */

#define BITFIRSTNZ(x) bitfirst(x)

#define BITFIRST(id, x) { dgword_t b_; BITFIRSTLOW(id, x, b_); }

/* index of nonzero bit */
INLINE int bitfirst(dgword_t x)
{
  dgword_t b;
  return bitfirstlow(x, &b);
}

#endif /* defined(__INTEL_COMPILER) && (DG_WORDBITS == 32) */

/* ********************** BITFIRST ENDS ****************************** */





typedef struct {
  int n;
  dgword_t *c; /* the jth bit of c[i] is 1 if the vertices i and j are adjacent */
} dg_t;



/* define a stock diagram object, with g->c declared by values `g_c' */
#define DG_MKSTOCK_(g, storeclass) \
  storeclass dgword_t g##_c[DG_NMAX]; \
  storeclass dg_t g[1] = {{ DG_NMAX, g##_c}};

#define DG_MKSTATICSTOCK(g) DG_MKSTOCK_(g, static)



/* open a diagram */
INLINE dg_t *dg_open(int n)
{
  dg_t *g;

#ifdef N
  die_if (n != N, "this version only support fixed N %d != n %d", N, n);
#endif
  xnew(g, 1);

  g->n = n;
  /* NOTE: some functions may not work for n == DG_NMAX */
  die_if (n > DG_NMAX, "do not support %d atoms\n", n);
  xnew(g->c, n);
  return g;
}



/* close the database */
INLINE void dg_close(dg_t *g)
{
  if (g->c) free(g->c);
  free(g);
}



/* copy diagram */
INLINE dg_t *dg_copy(dg_t *a, const dg_t *b)
{
  int i;
  DG_DEFN_(b)

  for (i = 0; i < DG_N_; i++)
    a->c[i] = b->c[i];
  return a;
}



INLINE dg_t *dg_clone(const dg_t *b)
{
  dg_t *a = dg_open(b->n);
  return dg_copy(a, b);
}



/* return if the edge between i and j are connected */
#define dg_linked(g, i, j) (int) ( ( (g)->c[i] >> (j) ) & 1u )



/* add an edge between i and j */
#define DG_LINK(g, i, j) { \
  (g)->c[i] |= MKBIT(j); \
  (g)->c[j] |= MKBIT(i); }

/* safer function version */
INLINE void dg_link(dg_t *g, int i, int j)
{
  DG_LINK(g, i, j);
}



/* remove an edge between i and j */
#define DG_UNLINK(g, i, j) { \
  (g)->c[i] &= ~MKBIT(j); \
  (g)->c[j] &= ~MKBIT(i); }

/* safer function version */
INLINE void dg_unlink(dg_t *g, int i, int j)
{
  DG_UNLINK(g, i, j);
}



/* remove all edges */
INLINE void dg_empty(dg_t *g)
{
  int i;
  DG_DEFN_(g)

  for (i = 0; i < DG_N_; i++)
    g->c[i] = 0;
}



/* a fully connected diagram */
INLINE void dg_full(dg_t *g)
{
  int i;
  DG_DEFN_(g)
  DG_DEFMASKN_()

  for (i = 0;  i < DG_N_; i++)
    /* all bits, except the ith, are 1s */
    g->c[i] = DG_MASKN_ ^ MKBIT(i);
}



/* the complement diagram of h: link <---> g: unlink */
INLINE dg_t *dg_complement(dg_t *h, dg_t *g)
{
  int i;
  DG_DEFN_(g)
  DG_DEFMASKN_();

  for (i = 0; i < DG_N_; i++)
    h->c[i] = DG_MASKN_ ^ (g->c[i] | MKBIT(i));
  return h;
}



#define dg_print0(g, nm)  dg_fprint0(g, stdout, nm)
#define dg_print(g)       dg_print0(g, NULL)

INLINE void dg_fprint0(const dg_t *g, FILE *fp, const char *nm)
{
  int i, j;
  DG_DEFN_(g)

  fprintf(fp, "%-4s", nm ? nm : "");
  for (i = 0; i < DG_N_; i++)
    fprintf(fp, " %2d", i);
  fprintf(fp, "\n");

  for (i = 0; i < DG_N_; i++) {
    fprintf(fp, "%2d: ", i);
    for (j = 0; j < DG_N_; j++)
      fprintf(fp, "  %c", dg_linked(g, i, j) ? '*' : ' ');
    fprintf(fp, "\n");
  }
}



/* degree of vertex i */
#define dg_deg(g, i) bitcount( (g)->c[i] )

/* degree of vertex i with respect to a vertex set vs */
#define dg_degvs(g, i, vs) bitcount( ((g)->c[i]) & vs )



/* get the degree sequence, return the number of edges */
INLINE int dg_degs(const dg_t *g, int *degs)
{
  int i, sum = 0;
  DG_DEFN_(g)

  for (i = 0; i < DG_N_; i++)
    sum += ( degs[i] = dg_deg(g, i) );
  return sum / 2;
}



/* count sort the degree sequence, TO BE TESTED */
INLINE void dg_csortdegs(int n, int *degs)
{
  int i, j, k = 0;
  static int cnt[DG_NMAX] = {0};
#pragma omp threadprivate(cnt)

  for (i = 0; i < DG_N_; i++) cnt[ degs[i] ]++;
  for (i = DG_N_ - 1; i >= 0; i--)
    for (j = 0; j < cnt[i]; j++) degs[k++] = i;
}



/* count the number of edges */
INLINE int dg_nedges(const dg_t *g)
{
  int i, ne;
  DG_DEFN_(g)

  for (ne = i = 0; i < DG_N_; i++)
    ne += bitcount(g->c[i]);
  return ne / 2; /* divided by 2 for double-counting */
}



/* check if a diagram is connected */
#define dg_connected(g) dg_connectedvs((g), MKBITSMASK(DG_GN_(g)))

/* check if the subgraph of the vertex set `vs' is connected */
INLINE int dg_connectedvs(const dg_t *g, dgword_t vs)
{
  dgword_t stack = vs & (-vs), bk; /* add the first vertex into the stack */
  int k;

  while (stack) {
    BITFIRSTLOW(k, stack, bk); /* first vertex (1-bit) in the stack */
    vs ^= bk; /* remove the vertex `k' from the to-do list `vs' */
    stack = (stack | g->c[k]) & vs; /* add neighbors of `k', limit to the to-do list */
    if (stack == vs) /* all rest vertices are connected to the tree */
      return 1;
  }
  return 0;
}



/* check if diagram is biconnected, bitwise version
 * we do not use lookup table by default */
INLINE int dg_biconnected(const dg_t *g)
{
  DG_DEFN_(g)
  DG_DEFMASKN_()
  dgword_t b;

#if !defined(N) || N <= 2
  if (DG_N_ > 2) {
#endif /* !defined(N) || N <= 2 */
    for (b = 1; b & DG_MASKN_; b <<= 1) {
      if ( !dg_connectedvs(g, b ^ DG_MASKN_) )
        return 0;
    }
    return 1;
#if !defined(N) || N <= 2
  } else if (DG_N_ == 2) {
    return (int) (g->c[1] & 1u);
  } else /* if (n < 2) */ {
    return 1;
  }
#endif /* !defined(N) || N <= 2 */
}



/* check if a sub-diagram is biconnected  */
INLINE int dg_biconnectedvs(const dg_t *g, dgword_t vs)
{
  dgword_t b, todo;

  for (todo = vs; todo; todo ^= b) {
    b = todo & (-todo); /* first vertex (1-bit) of the todo list */
    if ( !dg_connectedvs(g, vs ^ b) )
      return 0;
  }
  return 1;
}



/* check if a graph is biconnected
 * standard algorithm used as a reference
 * slower than the default bitwise version */
INLINE int dg_biconnected_std(const dg_t *g)
{
  int i0, v, par, id, root = 0;
  DG_DEFN_(g)
  static int stack[DG_NMAX + 1], parent[DG_NMAX], dfn[DG_NMAX], low[DG_NMAX];
#pragma omp threadprivate(stack, parent, dfn, low)

  for (v = 0; v < DG_N_; v++) {
    dfn[v] = 0; /* no vertex is visited */
    parent[v] = -1; /* no parent */
  }
  dfn[root] = low[root] = id = 1;
  stack[0] = root;
  stack[1] = -1;
  /* depth-first search to construct the spanning tree */
  for (i0 = 1; i0; ) {
    par = (i0 >= 1) ? stack[i0 - 1] : -1;
    if (i0 == DG_N_) {
      v = DG_N_; /* last level, no need to DFS */
    } else {
      /* v == the first unvisited vertex adjacent to `par' */
      for (v = ++stack[i0]; v < DG_N_; stack[i0] = ++v)
        /* dfn[] of an unvisited vertex is 0 */
        if ( dfn[v] == 0 && dg_linked(g, par, v) ) {
          dfn[v] = low[v] = ++id;
          parent[v] = par;
          stack[++i0] = -1; /* push `v', clear the next level */
          break; /* break the loop to go to the next level */
        }
    }
    if (v == DG_N_) { /* all children of `par' are visited, ready to pop */
      for (v = 0; v < DG_N_; v++) { /* update lowpoints */
        /* if we get here with i0 == 1, it means no vertex is connected to `root'
           if we get here with i0 == 2, it means there is vertex not connected
             to the first branch of the tree grown from `root', the graph is
             either unconnected or not biconnected */
        if ( i0 <= 2 && dfn[v] == 0 ) return 0;
        if ( !dg_linked(g, par, v) ) continue;
        if ( low[v] < low[par] && v != parent[par] )
          low[par] = low[v];
        /* the i0 == 1 case means `par' is the root, so it must be excluded */
        if ( i0 > 1 && low[v] >= dfn[par] && parent[v] == par )
          return 0;
      }
      i0--; /* pop */
    }
  }
  return 1;
}



/* code the connectivity */
INLINE dgword_t *dg_encode(const dg_t *g, dgword_t *code)
{
  int i, ib = 0;
  DG_DEFN_(g)
  dgword_t *c = code, ci;

  for (*c = 0, i = 0; i < DG_N_ - 1; i++) {
    *c |= ((ci = g->c[i]) >> (i + 1)) << ib;
    ib += DG_N_ - 1 - i;
#if !defined(N) || (N*(N-1)/2 > DG_WORDBITS)
    /* avoid the following step, if the number of edges
     * can be always contained in a single dgword_t */
    if (ib >= DG_WORDBITS) {
      ib -= DG_WORDBITS;
      *(++c) = ci >> (DG_N_ - ib);
    }
#endif
  }
  return code;
}



/* code the connectivity */
INLINE dg_t *dg_decode(dg_t *g, dgword_t *code)
{
  int i, j, ib = 0;
  DG_DEFN_(g)
  dgword_t *c = code;

  dg_empty(g);
  for (i = 0; i < DG_N_ - 1; i++) {
    for (j = i + 1; j < DG_N_; j++) {
      if ((*c >> ib) & 1u)
        DG_LINK(g, i, j);
      ++ib;
#if !defined(N) || (N*(N-1)/2 > DG_WORDBITS)
      /* avoid the following step, if the number of edges
       * can be always contained in a single dgword_t */
      if (ib == DG_WORDBITS) {
        ib = 0;
        c++;
      }
#endif
    }
  }
  return g;
}



/* compare two words of length 1 */
INLINE int dgword_cmp1(const dgword_t *a, const dgword_t *b)
{
  if (a[0] > b[0]) return 1;
  else if (a[0] < b[0]) return -1;
  return 0;
}



/* compare two words of length 2 */
INLINE int dgword_cmp2(const dgword_t *a, const dgword_t *b)
{
  if (a[0] > b[0]) return 1;
  else if (a[0] < b[0]) return -1;
  return 0;
}



/* compare two words of length n */
INLINE int dgword_cmp(const dgword_t *a, const dgword_t *b, int n)
{
  int k;
  for (k = n - 1; k >= 0; k--) {
    if (a[k] > b[k]) return 1;
    else if (a[k] < b[k]) return -1;
  }
  return 0;
}



/* check if two words are equal */
INLINE int dgword_eq1(const dgword_t *a, const dgword_t *b)
{
  return (a[0] == b[0]);
}



/* check if two 2-words are equal */
INLINE int dgword_eq2(const dgword_t *a, const dgword_t *b)
{
  return (a[0] == b[0]) && (a[1] == b[1]);
}



/* check if two n-words are equal */
INLINE int dgword_eq(const dgword_t *a, const dgword_t *b, int n)
{
  int k;
  for (k = 0; k < n; k++)
    if (a[k] != b[k]) return 0;
  return 1;
}



INLINE dgword_t *dgword_cpy(dgword_t *a, const dgword_t *b, int n)
{
  int k;

  for (k = 0; k < n; k++)
    a[k] = b[k];
  return a;
}



/* signed and unsigned comparison, and copy */
#if defined(N) && (N*(N-1)/2 <= DG_WORDBITS)
  #define dgcode_cmp(a, b, n) dgword_cmp1(a, b)
  #define dgcode_eq(a, b, n)  dgword_eq1(a, b)
  #define DGCODE_CPY(a, b, n) a[0] = b[0]
#elif defined(N) && (N*(N-1)/2 <= DG_WORDBITS * 2)
  #define dgcode_cmp(a, b, n) dgword_cmp2(a, b)
  #define dgcode_eq(a, b, n)  dgword_eq2(a, b)
  #define DGCODE_CPY(a, b, n) { a[0] = b[0]; a[1] = b[1]; }
#else
  #define dgcode_cmp(a, b, n) dgword_cmp(a, b, n)
  #define dgcode_eq(a, b, n)  dgword_eq(a, b, n)
  #define DGCODE_CPY(a, b, n) dgword_cpy(a, b, n)
#endif


#define dg_strset1(c) dg_sprintset(NULL, &(c), 1)
#define dg_sprintset1(s, c) dg_sprintset(s, &(c), 1)

INLINE char *dg_sprintset(char *s, const dgword_t *c, int nc)
{
  int ip, ic, k;
  char *p;
  dgword_t vs, b;

  /* allocate a string buffer if the input is empty */
  if (s == NULL && (s = malloc(nc * DG_WORDBITS * 20)) == NULL)
    return NULL;
  s[0] = '\0';
  p = s;
  for (ic = 0; ic < nc; ic++) {
    for (vs = c[ic]; vs; vs ^= b) {
      BITFIRSTLOW(k, vs, b);
      ip = sprintf(p, "%d ", DG_WORDBITS * ic + k);
      p += ip;
    }
  }
  if (p != s) /* remove the trailing space, if something has been written */
    p[-1] = '\0';
  return s;
}



#endif /* DG_H__ */

