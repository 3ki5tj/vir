#ifndef DG_H__
#define DG_H__
/* handling diagrams in the virial expansion
 * using bitwise operations
 * This file contains generic code only
 * one-word and multiple-word codes are
 * distributed in dg1.h and dgl.h */


#ifdef __INTEL_COMPILER
  /* 161: don't complain unknown omp pragma
   * 869: unreferenced parameter (n)
   * 981: operands evaluated in unspecified order
   * 2415: unused static variables */
  #pragma warning disable 161 869 981 2415
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



/* DG_WORDBITS used to be called CODEBITS */
#if !defined(DG_WORDBITS) && defined(CODEBITS)
#define DG_WORDBITS CODEBITS
#endif



#ifdef N
/* the user should not define NMAX and N */
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
#define DEFAULT_BITS DG_WORDBITS
#define ZCOM_PICK
#define ZCOM_UTIL
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



/* decide whether to use the one-word or multiple-word version */
#ifndef DGVS_ONEWORD
  #if DG_NMAX <= DG_WORDBITS
    #define DGVS_ONEWORD 1
  #else
    #define DGVS_ONEWORD 0
  #endif /* DG_NMAX <= DG_WORDBITS */
#endif /* !defined(DGVS_ONEWORD) */



#if DGVS_ONEWORD /* 1-word vertex set */

/* 1-word vertex set */
typedef dgword_t dgvs_t;
typedef dgword_t dgvsref_t;

#define DGVS_CLEAR(vs)              (vs) = 0;

#define dgvs_nonzero(vs)            ((vs) != 0)

/* vs = ~vs */
#define DGVS_NOT(vs)                (vs) = ~(vs);

/* copy vertex sets */
#define DGVS_CPY(vs1, vs2)          (vs1) = (vs2);

/* count the number of bits in the vertex set `vs' */
#define dgvs_count(vs)              bitcount(vs)

/* define a variable `iq' to be used in DGVS_FIRSTBIT1() or DGVS_FIRSTLOW1()
 * in the one-word case, nothing needs to be defined */
#define DGVS_DEFIQ_(iq)

/* save the bit corresponding to the first vertex to `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTBIT1(vs, b, iq)     (b) = (vs) & (-(vs));

/* save the index of the first vertex to `id'
 * and the bit corresponding to `id' is `b'
 * `vs' is dgvs_t, 'b' is dgword_t, `iq' is unused */
#define DGVS_FIRSTLOW1(k, vs, b, iq)  BITFIRSTLOW(k, vs, b)

#define dgvs_eq(a, b)           ((a) == (b))

/* a ^= b */
#define DGVS_XOR(a, b)          ((a) ^= (b));
#define DGVS_XOR1(a, b, iq)     DGVS_XOR(a, b)

/* a |= b */
#define DGVS_OR(a, b)           ((a) |= (b));
#define DGVS_OR1(a, b, iq)      DGVS_OR(a, b)

/* a &= b */
#define DGVS_AND(a, b)          ((a) &= (b));
#define DGVS_AND1(a, b, iq)     DGVS_AND(a, b)

/* a = b & c */
#define DGVS_AND2(a, b, c)      (a) = (b) & (c);

/* a &= ~b */
#define DGVS_ANDNOT(a, b)       ((a) &= ~(b));
#define DGVS_ANDNOT1(a, b, iq)  DGVS_AND(a, b)

/* add vertex `i' into the vertex set `vs' */
#define DGVS_ADD(vs, i)         DGVS_OR(vs, MKBIT(i))

/* remove vertex `i' from the vertex set `vs' */
#define DGVS_REMOVE(vs, i)      DGVS_ANDNOT(vs, MKBIT(i))

/* flip the state of vertex `i' in the vertex set `vs' */
#define DGVS_FLIP(vs, i)        DGVS_XOR(vs, MKBIT(i))

/* test if the vertex set `vs' has vertex `i' */
#define DGVS_HAS(vs, i)         (int) (((vs) >> i) & 0x1)

#define DGVS_MKBIT(i, b, iq)    b = MKBIT(i);

/* a mask with the lowest n words */
#define DGVS_MKBITSMASK(vs, n)  (vs) = MKBITSMASK(n);

#ifdef N
  #define DGVS_MASKN_ MKBITSMASK(N)
  #define DGVS_DEFMASKN_()
#else
  #define DGVS_MASKN_             maskn
  #define DGVS_DEFMASKN_()        dgvs_t maskn = MKBITSMASK(n);
#endif

/* a = ~b & set */
#define DGVS_COMPLM2(a, b, set) (a) = ~(b) & set;

#define dgvs_complement2(a, b, set) ((a) = (~(b) & set))



#else /* !DGVS_ONEWORD */



#if DG_WORDBITS == 32
  #define DG_NW         ((DG_NMAX + 31) >> 5)
  #define DG_IQ(i)      ((i) >> 5)
  #define DG_IR(i)      ((i) & 0x1f)
  #define DG_IB(i)      MKBIT( DG_IR(i) )
#elif DG_WORDBITS == 64
  #define DG_NW         ((DG_NMAX + 63) >> 6)
  #define DG_IQ(i)      ((i) >> 6)
  #define DG_IR(i)      ((i) & 0x3f)
  #define DG_IB(i)      MKBIT( DG_IR(i) )
#else
  #define DG_NW         ((DG_NMAX + DG_WORDBITS - 1) / DG_WORDBITS)
  #define DG_IQ(i)      ((i) / DG_WORDBITS)
  #define DG_IR(i)      ((i) % DG_WORDBITS)
  #define DG_IB(i)      MKBIT( DG_IR(i) )
#endif

/* multiple-word vertex set */
typedef dgword_t dgvs_t[DG_NW];
typedef dgword_t *dgvsref_t;

/* clear a vertex set */
#define DGVS_CLEAR(vs) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    (vs)[ii_] = 0; }

/* test if a vertex set is zero */
int dgvs_nonzero(dgvs_t vs)
{
  int i;
  for (i = 0; i < DG_NW; i++)
    if (vs[i] != 0) return 1;
  return 0;
}

/* negate a vertex set */
#define DGVS_NOT(vs) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    (vs)[ii_] = ~( (vs)[ii_] ); }

/* count the number of bits in the vertex set `vs' */
INLINE int dgvs_count(dgvs_t vs)
{
  int i, cnt = 0;

  for (i = 0; i < DG_NW; i++) cnt += bitcount(vs[i]);
  return cnt;
}

/* define a variable `iq' to be used in DGVS_FIRSTBIT1() or DGVS_FIRSTLOW1() */
#define DGVS_DEFIQ_(iq) int iq;

/* save the bit corresponding to the first vertex to `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTBIT1(vs, b, iq) { \
  dgword_t w_; \
  (b) = 0; \
  for (iq = 0; iq < DG_NW; iq++) { \
    if ((w_ = (vs)[iq]) != 0) { \
      (b) = w_ & (-w_); \
      break; \
    } } }

/* save the index of the first vertex to `id'
 * and the bit corresponding to `id' is `b'
 * `vs' is dgvs_t, `b' is dgword_t, `iq' is the word offset */
#define DGVS_FIRSTLOW1(k, vs, b, iq) { \
  dgword_t w_; \
  b = 0; k = 0; \
  for (iq = 0; iq < DG_NW; iq++) { \
    if ((w_ = (vs)[iq]) != 0) { \
      BITFIRSTLOW(k, w_, b); \
      k += iq * DG_WORDBITS; \
      break; \
    } } }

/* return if a == b */
INLINE int dgvs_eq(dgvs_t a, dgvs_t b)
{
  int i;
  for (i = 0; i < DG_NW; i++)
    if (a[i] != b[i]) return 0;
  return 1;
}

/* a ^= b */
#define DGVS_XOR(a, b) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) (a)[ii_] ^= (b)[ii_]; }

#define DGVS_XOR1(a, b, iq)     (a)[iq] ^= b;

/* a |= b */
#define DGVS_OR(a, b) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) (a)[ii_] |= (b)[ii_]; }

#define DGVS_OR1(a, b, iq)      (a)[iq] |= b;

/* a &= b */
#define DGVS_AND(a, b) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    (a)[ii_] &= (b)[ii_]; }

#define DGVS_AND1(a, b, iq)     (a)[iq] &= b;

/* a = b & c */
#define DGVS_AND2(a, b, c) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    (a)[ii_] = (b)[ii_] & (c)[ii_]; }

/* a &= ~b */
#define DGVS_ANDNOT(a, b) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    (a)[ii_] &= ~((b)[ii_]); }

#define DGVS_ANDNOT1(a, b, iq)  (a)[iq] &= ~b;

/* add vertex `i' into the vertex set `vs' */
#define DGVS_ADD(vs, i)         DGVS_OR1(vs, DG_IB(i), DG_IQ(i))

/* remove vertex `i' from the vertex set `vs' */
#define DGVS_REMOVE(vs, i)      DGVS_ANDNOT1(vs, DG_IB(i), DG_IQ(i))

/* flip the state of vertex `i' in the vertex set `vs' */
#define DGVS_FLIP(vs, i)        DGVS_XOR1(vs, DG_IB(i), DG_IQ(i))

/* test if the vertex set `vs' has vertex `i' */
#define DGVS_HAS(vs, i)         (int) (((vs)[DG_IQ(i)] >> DG_IR(i)) & 0x1)

/* copy a vertex set */
#define DGVS_CPY(vs1, vs2) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    vs1[ii_] = vs2[ii_]; }

#define DGVS_MKBIT(i, b, iq)    iq = DG_IQ(i), b = MKBIT(DG_IR(i));

/* a mask with the lowest n words */
#define DGVS_MKBITSMASK(vs, n) { int ii_; \
  int iq_ = DG_IQ(n), ir_ = DG_IR(n); \
  for (ii_ = 0; ii_ < iq_; ii_++) \
    (vs)[ii_] = MKINVBITSMASK(0); \
  if (ir_ != 0) (vs)[iq_++] = MKBITSMASK(ir_); \
  for (ii_ = iq_; ii_ < DG_NW; ii_++) \
    (vs)[ii_] = 0; }

#define DGVS_MASKN_       maskn
#define DGVS_DEFMASKN_()  dgvs_t maskn; DGVS_MKBITSMASK(maskn, DG_N_)

/* a = ~b & set */
#define DGVS_COMPLM2(a, b, set) { int ii_; \
  for (ii_ = 0; ii_ < DG_NW; ii_++) \
    a[ii_] = ~b[ii_] & set[ii_]; }

INLINE dgvsref_t dgvs_complement2(dgvs_t a, dgvs_t b, dgvs_t set)
{
  DGVS_COMPLM2(a, b, set)
  return a;
}

#endif /* DGVS_ONEWORD */

#define DG_MASKN_       DGVS_MASKN_
#define DG_DEFMASKN_    DGVS_DEFMASKN_




#define dgvs_print(vs, name) dgvs_fprint(vs, stdout, name)

/* print a vertex set */
INLINE void dgvs_fprint(dgvs_t vs, FILE *fp, const char *name)
{
  int i;

  if (name != NULL) fprintf(fp, "%-8s: [", name);
  for (i = 0; i < DG_NMAX; i++)
    fprintf(fp, "%c", DGVS_HAS(vs, i) ? '*' : ' ');
  fprintf(fp, "]\n");
}



#define dgvs_printn(vs, name) dgvs_fprintn(vs, stdout, name)

/* print a vertex set */
INLINE void dgvs_fprintn(dgvs_t vs, FILE *fp, const char *name)
{
  int v;
  dgword_t bv;
  dgvs_t c;
  DGVS_DEFIQ_(iq)

  if (name != NULL) fprintf(fp, "%-8s:", name);
  DGVS_CPY(c, vs)
  while ( dgvs_nonzero(c) ) {
    DGVS_FIRSTLOW1(v, c, bv, iq)
    fprintf(fp, " %d", v);
    DGVS_XOR1(c, bv, iq) /* remove vertex v */
  }
  fprintf(fp, "\n");
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

  DGVS_CPY(vs1, vs)
  DGVS_AND(vs1, g->c[i])
  return dgvs_count(vs1);
}

#endif /* DGVS_ONEWORD */



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
INLINE int dg_biconnected(const dg_t *g)
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



#else /* !DGVS_ONEWORD */



/* check if the subdiagram of the vertex set `vs' is connected */
INLINE int dg_connectedvs(const dg_t *g, dgvs_t vs0)
{
  dgvs_t vs, stack;
  dgword_t b;
  int iq, k;

  DGVS_CPY(vs, vs0)
  DGVS_CLEAR(stack)
  DGVS_FIRSTBIT1(vs, b, iq)
  stack[iq] = b;
  while ( dgvs_nonzero(stack) ) {
    DGVS_FIRSTLOW1(k, stack, b, iq) /* first vertex (1-bit) in the stack */
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
INLINE int dg_biconnected(const dg_t *g)
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
INLINE int dg_biconnectedvs(const dg_t *g, dgvs_t vs)
{
  int iq;
  dgvs_t todo;
  dgword_t bi;

  DGVS_CPY(todo, vs)
  while ( dgvs_nonzero(todo) ) {
    DGVS_FIRSTBIT1(todo, bi, iq) /* first vertex (1-bit) of the todo list */
    vs[iq] ^= bi;
    if ( !dg_connectedvs(g, vs) )
      return 0;
    vs[iq] ^= bi;
    todo[iq] ^= bi;
  }
  return 1;
}



#endif /* DGVS_ONEWORD */



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

  for (i = 0; i < DG_N_ - 1; i++) {
    for (j = i + 1; j < DG_N_; j++) {
      *c |= MKBIT(ib);
      if (++ib == DG_WORDBITS) {
        ib = 0;
        c++;
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

