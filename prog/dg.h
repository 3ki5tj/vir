#ifndef DIAGRAM_H__
#define DIAGRAM_H__
/* handling diagrams in the virial expansion
 * using bitwise operations */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RV3
#include "zcom.h"
#include <time.h>



#ifndef CODEBITS
/* Note:
 * #define CODEBITS (sizeof(code_t) * 8)
 * doesn't work, because the compiler cannot compute sizeof() in time
 * */
#define CODEBITS 32
#endif



#if CODEBITS == 32
typedef uint32_t code_t;
#else
typedef uint64_t code_t;
#endif


/* we only support a graph with at most CODEBITS vertices */
#define DG_NMAX CODEBITS



typedef struct {
  int n;
  code_t *c; /* if two particles are connected */
} dg_t;



/* count the number of 1 bits in x
 * http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable */
INLINE int bitcount(code_t x)
{
/* if the lowest two bits of x are zeros, and x has n nonzero bits,
 * then x, x+1, x+2, x+3 have respectivly n, n+1, n+1, n+2 nonzero bits */
#define B2(n) n,     n+1,     n+1,     n+2
#define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
  /* how many bits in 0..255 */
  static const unsigned char bits[256] = { B6(0), B6(1), B6(1), B6(2) };
  unsigned char *p = (unsigned char *) &x;

#if CODEBITS == 32
  return bits[p[0]] + bits[p[1]] + bits[p[2]] + bits[p[3]];
#elif CODEBITS == 64
  return bits[p[0]] + bits[p[1]] + bits[p[2]] + bits[p[3]]
       + bits[p[4]] + bits[p[5]] + bits[p[6]] + bits[p[7]];
#endif
}



/* invert the lowest k bits */
INLINE code_t bitinvert(code_t x, int k)
{
  return x ^ ~(~0u << k);
}



const int bruijn_index_[32] =
  { 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};

#define BRUIJNID32(b) bruijn_index_[((b) * 0x077CB531u) >> 27]
#define BRUIJNID64(b) (((b) & 0xFFFFFFFFu) \
    ? BRUIJNID32((b) & 0xFFFFFFFFu) : 32 + BRUIJNID32((b) >> 16))
#if CODEBITS == 32
  #define BRUIJNID(b) BRUIJNID32(b)
#elif CODEBITS == 64
  #define BRUIJNID(b) BRUIJNID64(b)
#endif



/* find the index of the lowest 1 bit
 * on return *b is the nonzero bit
 * http://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightMultLookup
 * ``Using de Bruijn Sequences to Index 1 in a Computer Word''
 * by Charles E. Leiserson, Harald Prokof, and Keith H. Randall.
 * http://supertech.csail.mit.edu/papers/debruijn.pdf */
INLINE int bitfirstlow(code_t x, code_t *b)
{
  (*b) = x & -x; /* such that only the lowest 1-bit survives */
  return BRUIJNID(*b);
}



/* index of nonzero bit */
INLINE int bitfirst(code_t x)
{
  code_t b;
  return bitfirstlow(x, &b);
}



/* return if the edge between i and j are connected */
#define dg_linked(g, i, j) ( ( (g)->c[i] >> (j) ) & 1u )



/* add an edge between i and j */
#define DG_LINK(g, i, j) { \
  (g)->c[i] |= 1u << (j); \
  (g)->c[j] |= 1u << (i); }

/* safer function version */
INLINE void dg_link(dg_t *g, int i, int j)
{
  DG_LINK(g, i, j);
}



/* remove an edge between i and j */
#define DG_UNLINK(g, i, j) { \
  (g)->c[i] &= ~(1u << (j)); \
  (g)->c[j] &= ~(1u << (i)); }

/* safer function version */
INLINE void dg_unlink(dg_t *g, int i, int j)
{
  DG_UNLINK(g, i, j);
}


/* construct `sg' by removing vertex `i0' from `g'
 * `sg' and `g' can be the same */
INLINE dg_t *dg_shrink1(dg_t *sg, const dg_t *g, int i0)
{
  int i, is = 0, n = g->n;
  code_t maskl, maskh;

  maskl = (1 << i0) - 1;
  maskh = ~maskl;
  for (i = 0; i < n; i++) {
    if (i0 == i) continue;
    sg->c[is++] = (g->c[i] & maskl) | ((g->c[i] >> 1) & maskh);
  }
  return sg;
}


/* degree of vertex i */
#define dg_deg(g, i) bitcount( (g)->c[i] )



/* get the degree sequence, return the number of edges */
INLINE int dg_degs(const dg_t *g, int *degs)
{
  int i, sum = 0, n = g->n;

  for (i = 0; i < n; i++)
    sum += ( degs[i] = dg_deg(g, i) );
  return sum / 2;
}



/* count sort the degree sequence, TO BE TESTED */
INLINE void dg_csortdegs(const dg_t *g, int *degs)
{
  int i, j, k = 0, n = g->n, cnt[DG_NMAX] = {0};

  for (i = 0; i < n; i++) cnt[ degs[i] ]++;
  for (i = n - 1; i >= 0; i--)
    for (j = 0; j < cnt[i]; j++) degs[k++] = i;
}



/* count the number of edges */
INLINE int dg_nedges(const dg_t *g)
{
  int i, ne;

  for (ne = i = 0; i < g->n; i++)
    ne += bitcount(g->c[i]);
  return ne / 2; /* divided by 2 for double-counting */
}



/* remove all edges */
INLINE void dg_empty(dg_t *g)
{
  int i, n = g->n;

  for (i = 0; i < n; i++)
    g->c[i] = 0;
}



/* a fully connected diagram */
INLINE void dg_full(dg_t *g)
{
  int i;
  code_t mask = (1 << g->n) - 1; /* the lowest n-bits are 1s */

  for (i = 0;  i < g->n; i++)
    /* all bits, except the ith, are 1s */
    g->c[i] = ~(1 << i) & mask;
}



/* the complement diagram of h: link <---> g: unlink */
INLINE dg_t *dg_complement(dg_t *h, dg_t *g)
{
  int i, n = g->n;
  code_t mask = (1 << g->n) - 1;

  for (i = 0; i < n; i++)
    h->c[i] = ~(g->c[i] | (1 << i)) & mask;
  return h;
}



/* code the connectivity */
INLINE code_t *dg_encode(const dg_t *g, code_t *code)
{
  int i, ib = 0, n = g->n;
  code_t *c = code, ci;

  for (*c = 0, i = 0; i < n - 1; i++) {
    *c |= ((ci = g->c[i]) >> (i + 1)) << ib;
    ib += n - 1 - i;
    if (ib >= CODEBITS) {
      ib -= CODEBITS;
      *(++c) = ci >> (n - ib);
    }
  }
  return code;
}



/* code the connectivity */
INLINE dg_t *dg_decode(dg_t *g, code_t *code)
{
  int i, j, ib = 0, n = g->n;
  code_t *c = code;

  dg_empty(g);
  for (i = 0; i < n - 1; i++) {
    for (j = i + 1; j < n; j++) {
      if ((*c >> ib) & 1u)
        DG_LINK(g, i, j);
      if (++ib == CODEBITS) {
        ib = 0;
        c++;
      }
    }
  }
  return g;
}



/* open a diagram */
static dg_t *dg_open(int n)
{
  dg_t *g;

  xnew(g, 1);
  g->n = n;
  die_if (n >= DG_NMAX, "do not support %d atoms\n", n);
  xnew(g->c, g->n);
  return g;
}



/* close the database */
static void dg_close(dg_t *g)
{
  if (g->c) free(g->c);
  free(g);
}



/* copy diagram */
INLINE dg_t *dg_copy(dg_t *a, const dg_t *b)
{
  int i;

  for (i = 0; i < b->n; i++)
    a->c[i] = b->c[i];
  return a;
}



INLINE dg_t *dg_clone(const dg_t *b)
{
  dg_t *a = dg_open(b->n);
  return dg_copy(a, b);
}



INLINE int dg_cmp(const dg_t *a, const dg_t *b)
{
  int i, df;

  for (i = 0; i < a->n; i++)
    if ((df = a->c[i] - b->c[i]) != 0)
      return df;
  return 0;
}



INLINE void dg_print(const dg_t *g)
{
  int i, j;

  printf("    ");
  for (i = 0; i < g->n; i++)
    printf(" %2d", i);
  printf("\n");

  for (i = 0; i < g->n; i++) {
    printf("%2d: ", i);
    for (j = 0; j < g->n; j++)
      printf("  %c", dg_linked(g, i, j) ? '*' : ' ');
    printf("\n");
  }
}


/* check if a diagram is connected */
#define dg_connected(g) dg_connectedvs(g, ((code_t) 1u << g->n) - 1u)

/* check if the subdiagram of the vertex set `vs' is connected */
INLINE int dg_connectedvs(const dg_t *g, code_t vs)
{
  code_t stack = vs & (-vs), b;

  while (stack) {
    int k = bitfirstlow(stack, &b); /* first vertex (1-bit) in the stack */
    vs ^= b; /* remove the vertex (1-bit) from the to-do list */
    stack = (stack | g->c[k]) & vs; /* update the stack */
    if ((stack ^ vs) == 0) return 1;
  }
  return 0;
}



/* check if diagram is biconnected, bitwise version
 * we do not use lookup table by default */
INLINE int dg_biconnected(const dg_t *g)
{
  int i, n = g->n;

  if (n > 2) {
    for (i = 0; i < n; i++) {
      code_t vs = (((code_t) 1u << n) - 1u) ^ ((code_t) 1u << i);
      if ( !dg_connectedvs(g, vs) )
        return 0;
    }
    return 1;
  } else if (n == 2) {
    return (g->c[0] >> 1) & 1u;
  } else /* if (n < 2) */ {
    return 1;
  }
}



/* check if a sub-diagram is biconnected  */
INLINE int dg_biconnectedvs(const dg_t *g, code_t vs)
{
  code_t b, todo;

  for (todo = vs; todo; todo ^= b) {
    b = todo & (-todo); /* first vertex (1-bit) of the todo list */
    if ( !dg_connectedvs(g, vs ^ b) )
      return 0;
  }
  return 1;
}



typedef short unqid_t;

typedef struct {
  int ng; /* number of unique diagrams */
  unqid_t *map; /* map a diagram to the unique diagram
                 `short' works for n <= 8 */
  code_t *first; /* index of the first unique diagram */
} dgmap_t;


/* static diagram map for n <= DGMAP_NMAX */
#define DGMAP_NMAX 8
dgmap_t dgmap_[DGMAP_NMAX + 1];



/* compute all permutations of n */
INLINE int dgmap_getperm(int n, int **pp)
{
  int np, npp, ipp, i;
  int top, st[DGMAP_NMAX + 2], used[DGMAP_NMAX + 2] = {0};

  for (np = 1, i = 1; i <= n; i++) np *= i;
  npp = np * n;
  xnew(*pp, npp);
  /* add all permutations of the diagram */
  for (i = 0; i < n; i++) {
    st[i] = -1;
    used[i] = 0;
  }
  for (ipp = 0, top = 0; ; ) {
    if (top >= n) {
      /* found an index permutation, build the graph */
      die_if (ipp >= npp, "ipp %d >= npp %d, n %d\n", ipp, npp, n);
      for (i = 0; i < n; i++) {
        die_if (st[i] < 0 || st[i] >= n, "bad %d/%d %d\n", i, n, st[i]);
        (*pp)[ipp++] = st[i];
      }
    } else {
      for (i = st[top]; ++i < n; )
        if ( !used[i] ) break;
      if (i < n) {
        used[i] = 1;
        st[top] = i;
        st[++top] = -1; /* clear the next level */
        continue;
      }
    }
    /* exhausted this level */
    if (--top < 0) break;
    used[ st[top] ] = 0;
  }
  return np;
}



/* compute initial diagram map */
INLINE int dgmap_init(dgmap_t *m, int n)
{
  code_t c, c1, ng, *masks, *ms;
  int ipm, npm, *pm, i, j, ipr, npr, sz;
  clock_t t0;

  die_if (n > DGMAP_NMAX, "n %d is too large\n", n);
  if (m->ng > 0) return 0; /* already initialized */

  t0 = clock();
  if (n >= 8) printf("n %d: initializing the diagram map\n", n);

  npr = n * (n - 1) / 2;
  ng = (code_t) 1u << npr;
  xnew(m->map, ng);
  for (c = 0; c < ng; c++) m->map[c] = -1;
  xnew(m->first, sz = 1024);

  /* compute all permutations */
  npm = dgmap_getperm(n, &pm);
  /* for each permutation compute the mask of each particle pair */
  xnew(masks, npm * npr);
  for (ipm = 0; ipm < npm; ipm++)
    for (ipr = 0, i = 0; i < n - 1; i++)
      for (j = i + 1; j < n; j++, ipr++)
        masks[ipm * npr + ipr] /* code bit of the pair (i, j) */
          = (code_t) 1u << getpairindex(pm[ipm*n + i], pm[ipm*n + j], n);
  free(pm);

  /* loop over all diagrams */
  for (m->ng = 0, c = 0; c < ng; c++) {
    if (m->map[c] >= 0) continue;
    if (m->ng >= sz) xrenew(m->first, sz += 1024);
    m->first[m->ng] = c;
    /* add all permutations of the diagram */
    for (ms = masks, ipm = 0; ipm < npm; ipm++) {
      /* `c1' is the code of the permutated diagram `c' */
      for (c1 = 0, ipr = 0; ipr < npr; ipr++, ms++)
        if ((c >> ipr) & 1u) c1 |= *ms;
      if (m->map[c1] < 0) m->map[c1] = m->ng;
      else die_if (m->map[c1] != m->ng,
        "error: corruption code: %#x %#x, graph %d, %d\n",
        c, c1, m->map[c1], m->ng);
    }
    m->ng++;
  }
  free(masks);
  printf("n %d, initialized, %d unique diagrams, %gs\n",
     n, m->ng, 1.*(clock() - t0)/CLOCKS_PER_SEC);
  return 0;
}



/* retrieve the diagram id */
INLINE unqid_t dg_getmapidx(const dg_t *g, code_t *c)
{
  int n = g->n;
  dgmap_t *m = dgmap_ + n;

  dgmap_init(m, n);
  dg_encode(g, c);
  return m->map[*c];
}



/*  retrieve the diagram id */
INLINE unqid_t dg_getmapid(const dg_t *g)
{
  code_t c;
  return dg_getmapidx(g, &c);
}



/* done */
INLINE void dgmap_done(void)
{
  int i;

  for (i = 1; i <= DGMAP_NMAX; i++) {
    dgmap_t *m = dgmap_ + i;
    if (m->map) free(m->map);
    if (m->first) free(m->first);
    m->ng = 0;
  }
}



/* the macro is slower than the direct version (due to dg_getmapid()) */
#define dg_biconnected_lookup(g) dg_biconnected_lookuplow(g->n, dg_getmapid(g))

/* check if a diagram is biconnected (lookup version)
 * use it only if `id' is known, otherwise it is slower */
INLINE int dg_biconnected_lookuplow(int n, unqid_t id)
{
  /* biconnectivity of unique diagrams */
  static int *bc[DGMAP_NMAX + 1] = {NULL};

  if (bc[n] == NULL) {
    dg_t *g1 = dg_open(n);
    dgmap_t *m = dgmap_ + n;
    int k;

    dgmap_init(m, n);
    xnew(bc[n], m->ng);
    for (k = 0; k < m->ng; k++) {
      code_t c = m->first[k];
      dg_decode(g1, &c);
      bc[n][k] = dg_biconnected(g1);
    }
    dg_close(g1);
  }
  return bc[ n ][ id ];
}



#endif
