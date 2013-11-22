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
#define NMAX N
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


#define DEFAULT_BITS DG_WORDBITS
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_BITS
#define ZCOM_RV3
#include "zcom.h"
#include <time.h>



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




#define mkbitsmask(n) MKBITSMASK(n) /* for compatibility */



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
    ne += dg_deg(g, i);
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
  #define dgcode_cmp(a, b, n) dgword_cmp1(a, b)
  #define dgcode_eq(a, b, n)  dgword_eq1(a, b)
  #define DG_CCPY(a, b, n) a[0] = b[0]
#elif defined(N) && (N*(N-1)/2 <= DG_WORDBITS * 2)
  #define dgcode_cmp(a, b, n) dgword_cmp2(a, b)
  #define dgcode_eq(a, b, n)  dgword_eq2(a, b)
  #define DG_CCPY(a, b, n) { a[0] = b[0]; a[1] = b[1]; }
#else
  #define dgcode_cmp(a, b, n) dgword_cmp(a, b, n)
  #define dgcode_eq(a, b, n)  dgword_eq(a, b, n)
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

