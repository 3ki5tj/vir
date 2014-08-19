#ifndef DG_H__
#define DG_H__
/* handling diagrams in the virial expansion
 * using bitwise operations */



#ifdef __INTEL_COMPILER
  /* 161: don't complain unknown omp pragma
   * 869: unreferenced parameter (n)
   * 981: operands evaluated in unspecified order
   * 2415: unused static variables */
  #pragma warning disable 161 869 981 2415
#elif defined(__GNUC__) && !defined(__INTEL_COMPILER)
  /* ignore all OpenMP pragmas */
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
  #pragma GCC diagonstic ignored "-Wunused-parameter"
#elif defined(_MSC_VER)
  /* 4068: unknown pragma
     4146: unary minus operator applied to unsigned type */
  #pragma warning(disable : 4068 4146)
  /* no need for "_s" versions of the CRT functions */
  #ifndef _CRT_SECURE_NO_DEPRECATE
  #define _CRT_SECURE_NO_DEPRECATE 1
  #endif
  #ifndef _CRT_SECURE_NO_WARNINGS
  #define _CRT_SECURE_NO_WARNINGS 1
  #endif
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



/* DG_WORDBITS used to be called CODEBITS */
#if !defined(DG_WORDBITS) && defined(CODEBITS)
#define DG_WORDBITS CODEBITS
#endif



#ifdef N
/* the user should not simultaneously define NMAX and N */
#define NMAX N
#endif /* defined(N) */



/* NMAX is a user shortcut for DG_NMAX
 * we shall use DG_NMAX below */
#if !defined(DG_NMAX) && defined(NMAX)
#define DG_NMAX NMAX
#endif



/* try to select a proper DG_WORDBITS if DG_NMAX or N is defined */
#if !defined(DG_WORDBITS) && defined(DG_NMAX)
  #if DG_NMAX <= 32
    #define DG_WORDBITS 32
  #elif DG_NMAX <= 64
    #define DG_WORDBITS 64
  #endif
#endif

/* set the default DG_WORDBITS */
#ifndef DG_WORDBITS
#define DG_WORDBITS 32
#endif /* defined(DG_WORDBITS) */




/* now that we have define DG_WORDBITS, we can use
 * the proper set of bitwise operations in zcom.h */
#define ZCOM_PICK
#define ZCOM_UTIL
#define WORDBITS_ DG_WORDBITS
#define ZCOM_BITS
#define ZCOM_RV3
#include "zcom.h"
#include <time.h>



/* define dgword_t and sdgword_t, these definitions must go
 * after we have included "zcom.h" to use uint32_t, etc */
#if DG_WORDBITS == 32
  typedef uint32_t  dgword_t;
  typedef int32_t   sdgword_t;
  #define DG_WPRI   PRIx32
  #define DG_WSCN   SCNx32
  #define DG_SWPRI  PRId32
  #define DG_SWSCN  SCNd32
#elif DG_WORDBITS == 64
  typedef uint64_t  dgword_t;
  typedef int64_t   sdgword_t;
  #define DG_WPRI   PRIx64
  #define DG_WSCN   SCNx64
  #define DG_SWPRI  PRId64
  #define DG_SWSCN  SCNd64
#else
  #error "bad DG_WORDBITS definition"
#endif




/* now that we have settle DG_WORDBITS, try to determine DG_NMAX
 * by default, we set DG_NMAX equal to DG_WORDBITS
 * so the 1-word version use is used */
#ifndef DG_NMAX
#define DG_NMAX DG_WORDBITS
#endif



/* if the user defined a fixed N, fixed length loops can be used
 * but all diagrams must be of size n
 * otherwise we define DG_N_ as n, whatever it is */
#ifdef  N
#define DG_N_ N
#define DG_GN_(g) N
#define DG_DEFN_(g)
#else
#define DG_N_ n
#define DG_GN_(g) ((g)->n)
#define DG_DEFN_(g) int n = (g)->n;
#endif
#define DG_NPR_ (DG_N_*(DG_N_ - 1)/2)



/* decide whether to use the one-word or multiple-word version */
#ifndef DGVS_ONEWORD
  #if DG_NMAX <= DG_WORDBITS
    #define DGVS_ONEWORD 1
  #else
    #define DGVS_ONEWORD 0
  #endif /* DG_NMAX <= DG_WORDBITS */
#endif /* !defined(DGVS_ONEWORD) */



#if DGVS_ONEWORD /* 1-word vertex set */

#define DG_NW 1

/* 1-word vertex set */
typedef dgword_t dgvs_t;
typedef dgword_t dgvsref_t;
#define DGVS_GETPTR(vs)   (&(vs)) /* dgvs_t --> (dgword_t *) */
#define DGVS_W2VS(w)      (w)     /* dgword_t --> dgvs_t */

#define DGVS_CLEAR(vs)              (vs) = 0;

#define dgvs_nonzero(vs)            ((vs) != 0)

#define dgvs_iszero(vs)             ((vs) == 0)

#define dgvs_neq(a, b)              ((a) != (b))

#define dgvs_eq(a, b)               ((a) == (b))

/* vs = ~vs */
#define DGVS_NOT(a, b)              (a) = (b) ^ MKBITSMASK(DG_N_);

/* copy vertex sets */
#define DGVS_CPY(a, b)              (a) = (b);

/* count the number of bits in the vertex set `vs' */
#define dgvs_count(vs)              bitcount(vs)

/* define a variable `iq' to be used in DGVS_FIRSTBIT() or DGVS_FIRSTLOW()
 * in the one-word case, nothing needs to be defined */
#define DGVS_DEFIQ_(iq)

#define DGVS_FIRSTWORD(vs)          (vs)

/* save the bit corresponding to the first vertex to `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTBIT(vs, b, iq)    (b) = (vs) & (-(vs));

/* save the index of the first vertex to `id'
 * and the bit corresponding to `id' is `b'
 * `vs' is dgvs_t, 'b' is dgword_t, `iq' is unused */
#define DGVS_FIRSTLOW(k, vs, b, iq)   BITFIRSTLOW(k, vs, b)

#define dgvs_first(vs)          bitfirst(vs)

#define dgvs_bit2id(vs)         BIT2ID(vs)

/* a ^= b */
#define DGVS_XOR(a, b)          ((a) ^= (b));
#define DGVS_XOR1(a, b, iq)     DGVS_XOR(a, b)

/* a = b ^ c */
#define DGVS_XOR2(a, b, c)      (a) = (b) ^ (c);
#define dgvs_xor2(a, b, c)      ((a) = (b) ^ (c))

/* a |= b */
#define DGVS_OR(a, b)           ((a) |= (b));
#define DGVS_OR1(a, b, iq)      DGVS_OR(a, b)

/* a = b | c */
#define DGVS_OR2(a, b, c)       (a) = (b) | (c);
#define dgvs_or2(a, b, c)       ((a) = (b) | (c))

/* a &= b */
#define DGVS_AND(a, b)          ((a) &= (b));
#define DGVS_AND1(a, b, iq)     DGVS_AND(a, b)

/* a = b & c */
#define DGVS_AND2(a, b, c)      (a) = (b) & (c);
#define dgvs_and2(a, b, c)      ((a) = (b) & (c))

/* a &= ~b */
#define DGVS_MINUS(a, b)        ((a) &= ~(b));
#define DGVS_MINUS1(a, b, iq)   DGVS_MINUS(a, b)

/* a = b & ~c */
#define DGVS_MINUS2(a, b, c)    (a) = (b) & ~(c);
#define dgvs_minus2(a, b, c)    ((a) = (b) & ~(c))

/* add vertex `i' into the vertex set `vs' */
#define DGVS_ADD(vs, i)         DGVS_OR(vs, MKBIT(i))

/* remove vertex `i' from the vertex set `vs' */
#define DGVS_REMOVE(vs, i)      DGVS_MINUS(vs, MKBIT(i))

/* flip the state of vertex `i' in the vertex set `vs' */
#define DGVS_FLIP(vs, i)        DGVS_XOR(vs, MKBIT(i))

/* test if the vertex set `vs' has vertex `i' */
#define DGVS_HAS(vs, i)         (int) (((vs) >> i) & 0x1)

/* right shift one bit */
#define DGVS_RSHIFT1(a, b)      (a) = (b) >> 1;

/* left shift one bit */
#define DGVS_LSHIFT1(a, b)      (a) = (b) << 1;

/* test if the vertex set `vs' has bit `bi' */
#define DGVS_HASBIT(vs, bi, iq) (((vs) & (bi)) != 0)

/* set `vs' as a vertex set of a single bit `bi' */
#define DGVS_ONEBIT(vs, bi, iq) (vs) = (bi);

#define DGVS_MKBIT(i, bi, iq)    bi = MKBIT(i);

/* a mask with the lowest n words */
#define DGVS_MKBITSMASK(vs, n)  (vs) = MKBITSMASK(n);

#define dgvs_mkbitsmask(vs, n)  ((vs) = MKBITSMASK(n))

#ifdef N
  #define DGVS_MASKN_             MKBITSMASK(N)
  #define DGVS_DEFMASKN_()
#else
  #define DGVS_MASKN_             maskn
  #define DGVS_DEFMASKN_()        dgvs_t maskn = MKBITSMASK(n);
#endif

#define DGVS_MKINVSET(vs, n, i) { vs = MKBITSMASK(n) ^ MKBIT(i); }
#define dgvs_mkinvset(vs, n, i) (vs = MKBITSMASK(n) ^ MKBIT(i))

#define DGVS_MKINVSET2(vs, n, i, j)  { vs = MKBITSMASK(n) ^ MKBIT(i) ^ MKBIT(j); }
#define dgvs_mkinvset2(vs, n, i, j)  (vs = MKBITSMASK(n) ^ MKBIT(i) ^ MKBIT(j))

#else /* !DGVS_ONEWORD */



#define DG_NW       BITS_GETNW(DG_NMAX)
#define DG_IQ(i)    BITS_GETNQ(i)
#define DG_IR(i)    BITS_GETNR(i)
#define DG_IB(i)    BITS_GETNB(i)

/* multiple-word vertex set */
typedef dgword_t dgvs_t[DG_NW];
typedef dgword_t *dgvsref_t;
#define DGVS_GETPTR(vs) (vs)    /* dgvs_t --> (dgword_t *) */
#define DGVS_W2VS(w)    (&(w))  /* dgword_t --> dgvs_t */

/* define a variable `iq' to be used in DGVS_FIRSTBIT() or DGVS_FIRSTLOW() */
#define DGVS_DEFIQ_(iq) int iq;

/* copy a vertex set */
#define DGVS_CPY(a, b) BITS_CPY(a, b, DG_NW)

/* clear a vertex set */
#define DGVS_CLEAR(vs) BITS_CLEAR(vs, DG_NW)

/* test if a vertex set is zero */
#define dgvs_nonzero(vs)  bits_nonzero(vs, DG_NW)
#define dgvs_iszero(vs)   bits_iszero(vs, DG_NW)

/* test if two vertex sets are equal */
#define dgvs_neq(a, b)    bits_notequal(a, b, DG_NW)
#define dgvs_eq(a, b)     bits_equal(a, b, DG_NW)

/* count the number of bits in the vertex set `vs' */
#define dgvs_count(vs)    bits_count(vs, DG_NW)

/* negate a vertex set */
#define DGVS_NOT(a, b)    BITS_NINV(a, b, DG_NMAX, DG_NW)

#define DGVS_FIRSTWORD(vs)  (vs)[0]

/* save the bit corresponding to the first vertex to `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTBIT(vs, b, iq) BITS_FIRSTBIT(vs, DG_NW, b, iq)

/* save the index of the first vertex to `id'
 * and the bit corresponding to `id' is `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTLOW(id, vs, b, iq) BITS_FIRSTLOW(id, vs, DG_NW, b, iq)

/* first vertex in `vs' */
#define dgvs_first(vs) bits_first(vs, DG_NW)

/* assuming `vs' contains a single bit, return the bit `id' */
#define dgvs_bit2id(vs) bits_bit2id(vs, DG_NW)

/* a ^= b */
#define DGVS_XOR(a, b)          BITS_XOR(a, b, DG_NW)
#define DGVS_XOR1(a, b, iq)     BITS_XOR1(a, b, iq)
/* a = b ^ c */
#define DGVS_XOR2(a, b, c)      BITS_XOR2(a, b, c, DG_NW)
#define dgvs_xor2(a, b, c)      bits_xor2(a, b, c, DG_NW)

/* a |= b */
#define DGVS_OR(a, b)           BITS_OR(a, b, DG_NW)
#define DGVS_OR1(a, b, iq)      BITS_OR1(a, b, iq)
/* a = b | c */
#define DGVS_OR2(a, b, c)       BITS_OR2(a, b, c, DG_NW)
#define dgvs_or2(a, b, c)       bits_or2(a, b, c, DG_NW)

/* a &= b */
#define DGVS_AND(a, b)          BITS_AND(a, b, DG_NW)
#define DGVS_AND1(a, b, iq)     BITS_AND1(a, b, iq)
/* a = b & c */
#define DGVS_AND2(a, b, c)      BITS_AND2(a, b, c, DG_NW)
#define dgvs_and2(a, b, c)      bits_and2(a, b, c, DG_NW)

/* a &= ~b, or a -= b (set complement) */
#define DGVS_MINUS(a, b)        BITS_MINUS(a, b, DG_NW)
#define DGVS_MINUS1(a, b, iq)   BITS_MINUS1(a, b, iq)
/* a = b & ~c or a = b - c (set complement) */
#define DGVS_MINUS2(a, b, c)    BITS_MINUS2(a, b, c, DG_NW)
#define dgvs_minus2(a, b, c)    bits_minus2(a, b, c, DG_NW)

/* add vertex `i' into the vertex set `vs' */
#define DGVS_ADD(vs, i)         BITS_ADD(vs, i)

/* remove vertex `i' from the vertex set `vs' */
#define DGVS_REMOVE(vs, i)      BITS_REMOVE(vs, i)

/* flip the state of vertex `i' in the vertex set `vs' */
#define DGVS_FLIP(vs, i)        BITS_FLIP(vs, i)

/* test if the vertex set `vs' has vertex `i' */
#define DGVS_HAS(vs, i)         BITS_HAS(vs, i)

/* test if the vertex set `vs' has bit `bi' at offset `iq' */
#define DGVS_HASBIT(vs, bi, iq) BITS_HASBIT(vs, bi, iq)

/* right shift one bit */
#define DGVS_RSHIFT1(a, b) BITS_RSHIFT1(a, b, DG_NW)

/* left shift one bit */
#define DGVS_LSHIFT1(a, b) BITS_LSHIFT1(a, b, DG_NW)

#define DGVS_MKBIT(i, bi, iq)    BITS_MKBIT(i, bi, iq)
/* set `vs' as a vertex set of a single bit `bi' */
#define DGVS_ONEBIT(vs, bi, iq)  BITS_ONEBIT(vs, DG_NW, bi, iq)

/* make a mask for the lowest k elements */
#define DGVS_MKBITSMASK(vs, k) BITS_MKBITSMASK(vs, k, DG_NW)
#define dgvs_mkbitsmask(vs, k) bits_mkbitsmask(vs, k, DG_NW)



#define DGVS_MASKN_       maskn
#define DGVS_DEFMASKN_()  dgvs_t maskn; DGVS_MKBITSMASK(maskn, DG_N_)

#define DGVS_MKINVSET(vs, n, i) BITS_MKINVSET(vs, n, DG_NW, i)
#define dgvs_mkinvset(vs, n, i) bits_mkinvset(vs, n, DG_NW, i)

#define DGVS_MKINVSET2(vs, n, i, j)  BITS_MKINVSET2(vs, n, DG_NW, i, j)
#define dgvs_mkinvset2(vs, n, i, j)  bits_mkinvset2(vs, n, DG_NW, i, j)

#endif /* DGVS_ONEWORD */



#define DG_MASKN_       DGVS_MASKN_
#define DG_DEFMASKN_    DGVS_DEFMASKN_



#define dgvs_print(vs, name) dgvs_fprint(vs, stdout, name)

#define dgvs_fprint(vs, fp, name) \
  dgvs_fprint0(vs, DG_NMAX, fp, name, "\n")

/* print a set as star or blank pattern */
INLINE void dgvs_fprint0(dgvs_t vs, int n, FILE *fp,
    const char *name, const char *ending)
{
  int i;

  if (fp == NULL) fp = stdout;
  if (name != NULL) fprintf(fp, "%-8s: ", name);
  fprintf(fp, "[");
  for (i = 0; i < n; i++)
    fprintf(fp, "%c", DGVS_HAS(vs, i) ? '*' : ' ');
  fprintf(fp, "]");
  if (ending) fputs(ending, fp);
}



#define dgvs_printn(vs, name) dgvs_fprintn(vs, stdout, name)

#define dgvs_fprintn(vs, fp, name) \
  dgvs_fprintn0(vs, DG_NMAX, fp, name, "\n")

/* print a set as numbers */
INLINE void dgvs_fprintn0(dgvs_t vs, int n, FILE *fp,
    const char *name, const char *ending)
{
  int i;

  if (fp == NULL) fp = stdout;
  if (name != NULL) fprintf(fp, "%-8s:", name);
  for (i = 0; i < n; i++)
    if ( DGVS_HAS(vs, i) ) fprintf(fp, " %d", i);
  if (ending) fputs(ending, fp);
}





typedef struct {
  int n;
  dgvs_t *c;
} dg_t;



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

  for (i = 0; i < DG_N_; i++) {
    DGVS_CPY(a->c[i], b->c[i])
  }
  return a;
}



/* compare diagrams */
INLINE int dg_cmp(const dg_t *a, const dg_t *b)
{
  int i;
  DG_DEFN_(b)

  for (i = 0; i < DG_N_; i++) {
    if (a->c[i] != b->c[i]) return a->c[i] > b->c[i] ? 1 : -1;
  }
  return 0;
}



/* remove all edges */
INLINE void dg_empty(dg_t *g)
{
  int i;
  DG_DEFN_(g)

  for (i = 0; i < DG_N_; i++)
    DGVS_CLEAR(g->c[i])
}



/* a fully connected diagram */
INLINE void dg_full(dg_t *g)
{
  int i;
  dgvs_t vs;
  DG_DEFN_(g)
  DGVS_DEFMASKN_()

  dg_empty(g);
  for (i = 0;  i < DG_N_; i++) {
    DGVS_CPY(vs, DGVS_MASKN_)
    DGVS_REMOVE(vs, i)
    DGVS_CPY(g->c[i], vs)
  }
}



#define dg_linked(g, i, j) DGVS_HAS((g)->c[i], j)


/* add an edge between i and j */
#define DG_LINK(g, i, j) { \
  DGVS_ADD((g)->c[i], j) \
  DGVS_ADD((g)->c[j], i) }

/* safer function version */
INLINE void dg_link(dg_t *g, int i, int j)
{
  DG_LINK(g, i, j);
}



/* remove an edge between i and j */
#define DG_UNLINK(g, i, j) { \
  DGVS_REMOVE((g)->c[i], j) \
  DGVS_REMOVE((g)->c[j], i) }

/* safer function version */
INLINE void dg_unlink(dg_t *g, int i, int j)
{
  DG_UNLINK(g, i, j);
}



/* degree of vertex i */
#define dg_deg(g, i) dgvs_count( (g)->c[i] )


/* degree of vertex i with respect to a vertex set vs */
#if DGVS_ONEWORD

#define dg_degvs(g, i, vs) dgvs_count( ((g)->c[i]) & vs )

#else /* !DGVS_ONEWORD */

INLINE int dg_degvs(const dg_t *g, int i, dgvs_t vs)
{
  dgvs_t vs1;
  return dgvs_count( dgvs_and2(vs1, g->c[i], vs) );
}

#endif /* DGVS_ONEWORD */



/* static variables for the degree sequence */
static int dg_nedges_;
static int dg_degs_[DG_NMAX];
#pragma omp threadprivate(dg_nedges_, dg_degs_)

/* compute the degrees of all vertices if needed
 * `ned' is (int *), if it is NULL on input, we point it to `&ned_local'
 * `degs' is (int *), if it is NULL on input, we point it to `degs_local'
 * if `ned' == NULL or if `*ned' <= 0, compute the degree sequence
 * otherwise, we assume that the values in `*ned' and `degs' are valid
 * and the two local variables are unused.
 * An example call:
 *  DG_CALC_DEGS(ned, degs, dg_nedges_, dg_degs_)
 * */
#define DG_CALC_DEGS(ned, degs, ned_local, degs_local) { \
  if ( ned == NULL ) { ned = &ned_local, ned_local = -1; } \
  if ( degs == NULL ) degs = degs_local; \
  if ( *ned <= 0 ) *ned = dg_degs(g, degs); \
}



/* get the degree sequence, return the number of edges */
INLINE int dg_degs(const dg_t *g, int *degs)
{
  int i, sum = 0;
  DG_DEFN_(g)

  for (i = 0; i < DG_N_; i++)
    sum += ( degs[i] = dg_deg(g, i) );
  return sum / 2;
}



/* count the number of edges */
INLINE int dg_nedges(const dg_t *g)
{
  int i, ne;
  DG_DEFN_(g)

  for (ne = i = 0; i < DG_N_; i++)
    ne += dg_deg(g, i);
  return ne / 2; /* divided by 2 for double counting */
}



#define dg_print0(g, nm)  dg_fprint0(g, stdout, nm)
#define dg_print(g)       dg_print0(g, NULL)
#define dg_fprint(g, fp)  dg_fprint0(g, fp, NULL)

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



#if DGVS_ONEWORD



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
INLINE int dg_biconnected_simple(const dg_t *g)
{
  DG_DEFN_(g)
  DGVS_DEFMASKN_()
  dgword_t b;

#if !defined(N) || N <= 2
  if (DG_N_ > 2) {
#endif /* !defined(N) || N <= 2 */
    for (b = 1; b & DGVS_MASKN_; b <<= 1) {
      if ( !dg_connectedvs(g, b ^ DGVS_MASKN_) )
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
INLINE int dg_biconnectedvs_simple(const dg_t *g, dgword_t vs)
{
  dgword_t b, todo;

  for (todo = vs; todo; todo ^= b) {
    b = todo & (-todo); /* first vertex (1-bit) of the to-do list */
    if ( !dg_connectedvs(g, vs ^ b) )
      return 0;
  }
  return 1;
}



#else /* !DGVS_ONEWORD */



/* check if the subdiagram of the vertex set `vs' is connected */
static int dg_connectedvs(const dg_t *g, dgvs_t vs0)
{
  dgvs_t vs, stack;
  dgword_t b;
  int iq, k;

  DGVS_CPY(vs, vs0)
  DGVS_FIRSTBIT(vs, b, iq)
  DGVS_CLEAR(stack)
  stack[iq] = b;
  while ( dgvs_nonzero(stack) ) {
    DGVS_FIRSTLOW(k, stack, b, iq) /* first vertex (1-bit) in the stack */
    DGVS_XOR1(vs, b, iq) /* remove the vertex (1-bit) from the to-do list */
    /* update the stack */
    DGVS_OR(stack, g->c[k]) /* stack |= q->c[k] */
    DGVS_AND(stack, vs) /* stack &= vs */
    /* if all vertices left are all in the stack, it is connected */
    if ( dgvs_eq(stack, vs) ) return 1;
  }
  return 0;
}



/* check if a diagram is connected */
INLINE int dg_connected(const dg_t *g)
{
  DG_DEFN_(g)
  DGVS_DEFMASKN_()
  return dg_connectedvs(g, maskn);
}



/* check if diagram is biconnected, bitwise version
 * we do not use lookup table by default */
INLINE int dg_biconnected_simple(const dg_t *g)
{
  int i, iq;
  DG_DEFN_(g)
  dgvs_t vs;
  dgword_t b;

  DGVS_MKBITSMASK(vs, DG_N_)
  for (i = 0; i < DG_N_; i++) {
    iq = DG_IQ(i);
    b = MKBIT( DG_IR(i) );
    vs[iq] ^= b;
    if ( !dg_connectedvs(g, vs) )
      return 0;
    vs[iq] ^= b;
  }
  return 1;
}



/* check if a sub-diagram is biconnected  */
INLINE int dg_biconnectedvs_simple(const dg_t *g, dgvs_t vs)
{
  int iq;
  dgvs_t todo;
  dgword_t bi;

  DGVS_CPY(todo, vs)
  while ( dgvs_nonzero(todo) ) {
    DGVS_FIRSTBIT(todo, bi, iq) /* first vertex (1-bit) of the to-do list */
    vs[iq] ^= bi;
    if ( !dg_connectedvs(g, vs) )
      return 0;
    vs[iq] ^= bi;
    todo[iq] ^= bi;
  }
  return 1;
}



#endif /* DGVS_ONEWORD */



#ifndef DGBC_SIMPLE
#define dg_biconnected(g)         dg_biconnected_std(g, 0)
#else
#define dg_biconnected(g)         dg_biconnected_simple(g)
#endif
#define dg_biconnectedvs(g, vs)   dg_biconnectedvs_simple(g, vs)



/* check if a graph is biconnected
 * standard algorithm with bitwise operation */
static int dg_biconnected_std(const dg_t *g, int root)
{
  int top, v, par, ppar, id;
  DG_DEFN_(g)
  static int stack[DG_NMAX + 1], dfn[DG_NMAX], low[DG_NMAX];
#pragma omp threadprivate(stack, dfn, low)
  dgvs_t unvisited, gc;
  dgword_t bv;
  DGVS_DEFIQ_(vq)

  dfn[root] = low[root] = id = 1;
  DGVS_MKINVSET(unvisited, DG_N_, root) /* all vertices except `root' */
  stack[0] = root;
  /* depth-first search to construct the spanning tree */
  for (top = 1; top != 0; ) {
    par = stack[top - 1];

    /* find the first unvisited vertex `v' adjacent to `par'
     * the set of the above vertices is `gc' */
    if ( dgvs_nonzero( dgvs_and2(gc, g->c[par], unvisited) ) ) {
      /* such a vertex exists, push */
      DGVS_FIRSTLOW(v, gc, bv, vq)
      DGVS_XOR1(unvisited, bv, vq) /* remove `bv' from `unvisited' */
      stack[top++] = v;
      dfn[v] = low[v] = ++id;
      continue; /* push */
    }
    /* if no vertex for this level, fall through and pop */

    /* if we get here with top == 1, we have finished enumerating all vertices
     *   connected to the root, if there is any unvisited vertex, it is not
     *   connected to `root' */
    /* if we get here with top == 2, we have finished enumerating all vertices
     *   int the first branch of the tree grown from `root', if there is any
     *   unvisited vertex, the graph is disconnected or not biconnected,
     *   in the latter case, `root' is an articulation point. */
    if (top <= 2) return dgvs_iszero(unvisited);

    /* all children of `par' are visited, ready to pop */
    /* update the low-point of `par', low[par] */
    DGVS_CPY(gc, g->c[par])
    ppar = stack[top - 2];
    DGVS_REMOVE(gc, ppar)
    while ( dgvs_nonzero(gc) ) {
      DGVS_FIRSTLOW(v, gc, bv, vq)
      DGVS_XOR1(gc, bv, vq)
      /* low[par] = min{low[par], low[v]}
       * if `v' is a son of `par' in the tree, this formula includes
       *    backedges that start from `v' or a descendant of `v'
       * otherwise, `v' holds a smaller `dfn' than `par', so
       *    par-->v is a backedge  */
      if ( low[v] < low[par] )
        low[par] = low[v];
    }
    //printf("popping top %d, par %d, low %d, ppar %d, stack %d, %d, %d\n",
    //    top, par, low[par], ppar, stack[0], stack[1], stack[2]); getchar();
    /* now low[par] is correctly computed */
    if ( low[par] >= dfn[ppar] )
      return 0;
    top--; /* pop */
  }
  return 1; /* should never reach here */
}



#if DGVS_ONEWORD



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



#else /* !DGVS_ONEWORD */



/* encode the connectivity */
INLINE dgword_t *dg_encode(const dg_t *g, dgword_t *code)
{
  int i, j, ib = 0;
  dgword_t *c = code;
  DG_DEFN_(g)

  *c = 0;
  for (i = 0; i < DG_N_ - 1; i++) {
    for (j = i + 1; j < DG_N_; j++) {
      if ( dg_linked(g, i, j) )
        *c |= MKBIT(ib);
      if (++ib == DG_WORDBITS) {
        ib = 0;
        *(++c) = 0;
      }
    }
  }
  return code;
}



/* code the connectivity */
INLINE dg_t *dg_decode(dg_t *g, dgword_t *code)
{
  int i, j, ib = 0;
  dgword_t *c = code;
  DG_DEFN_(g)

  dg_empty(g);
  for (i = 0; i < DG_N_ - 1; i++) {
    for (j = i + 1; j < DG_N_; j++) {
      if ((*c >> ib) & 1u)
        DG_LINK(g, i, j);
      if (++ib == DG_WORDBITS) {
        ib = 0;
        c++;
      }
    }
  }
  return g;
}



#endif /* DGVS_ONEWORD */



#define dg_printcode(c, n) dg_fprintcode(stdout, c, n)

/* print out the code */
INLINE void dg_fprintcode(FILE *fp, const dgword_t *c, int n)
{
  int i;
  for (i = 0; i < n; i++)
    fprintf(fp, "%#" DG_WPRI " ", c[i]);
  fprintf(fp, "\n");
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
  if (a[1] > b[1]) return 1;
  else if (a[1] < b[1]) return -1;
  if (a[0] > b[0]) return 1;
  else if (a[0] < b[0]) return -1;
  return 0;
}



/* compare two words of length 3 */
INLINE int dgword_cmp3(const dgword_t *a, const dgword_t *b)
{
  if (a[2] > b[2]) return 1;
  else if (a[2] < b[2]) return -1;
  if (a[1] > b[1]) return 1;
  else if (a[1] < b[1]) return -1;
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



/* the number of words to save the code of a graph of vertices
 * the number of bits needed is n * (n - 1) / 2
 * divide this by sizeof(dgword_t) is DG_CWORDS */
#define DG_CBITS  (DG_NMAX *(DG_NMAX - 1)/2)
#define DG_CWORDS ((DG_CBITS + DG_WORDBITS - 1)/DG_WORDBITS)
#define DG_CPRI   DG_WPRI
#define DG_CSCN   DG_WSCN


/* map DG_CWORDS_ to the macro if possible, or a variable otherwise */
#ifdef N
#define DG_CWORDS_(cwords) DG_CWORDS
#else
#define DG_CWORDS_(cwords) cwords
#endif



/* copy the code */
/* signed and unsigned comparison of codes, and copy */
#if defined(N) && (N*(N-1)/2 <= DG_WORDBITS)
  #define DG_CCMP(a, b, n) dgword_cmp1(a, b)
  #define DG_CEQ(a, b, n)  ((a)[0] == (b)[0])
  #define DG_CCPY(a, b, n) a[0] = b[0]
#elif defined(N) && (N*(N-1)/2 <= DG_WORDBITS * 2)
  #define DG_CCMP(a, b, n) dgword_cmp2(a, b)
  #define DG_CEQ(a, b, n)  ((a)[0] == (b)[0] && (a)[1] == (b)[1])
  #define DG_CCPY(a, b, n) { a[0] = b[0]; a[1] = b[1]; }
#elif defined(N) && (N*(N-1)/2 <= DG_WORDBITS * 3)
  #define DG_CCMP(a, b, n) dgword_cmp3(a, b)
  #define DG_CEQ(a, b, n)  \
    ((a)[0] == (b)[0] && (a)[1] == (b)[1] && (a)[2] == (b)[2])
  #define DG_CCPY(a, b, n) { a[0] = b[0]; a[1] = b[1]; a[2] = b[2]; }
#else
  #define DG_CCMP(a, b, n) dgword_cmp(a, b, n)
  #define DG_CEQ(a, b, n)  dgword_eq(a, b, n)
  #define DG_CCPY(a, b, n) dgword_cpy(a, b, n)
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

