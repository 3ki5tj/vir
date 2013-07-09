#ifndef DIAGRAM_H__
#define DIAGRAM_H__
/* handling diagrams in the virial expansion
 * using bitwise operations */


#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RV3
#include "zcom.h"

/* currently only support 32-bit */
typedef uint32_t code_t;

/* Note: 
 * #define CODEBITS (sizeof(code_t) * 8)
 * doesn't work
 * */
#define CODEBITS 32


typedef struct {
  int n;
  code_t *c; /* if two particles are connected */
} dg_t;



/* count the number of 1 bits in x */
INLINE int bitcount(code_t x)
{
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



/* find the lowest 1 bit */
INLINE int bitfirst(code_t x)
{
#if CODEBITS == 32
  static const int index[32] =
  { 0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};

  x &= -x; /* such that only the lowest 1-bit survives */
  return index[(code_t)(x * 0x077CB531u) >> 27];
#endif
}



/* return if the edge between i and j are connected */
INLINE int dg_linked(const dg_t *g, int i, int j)
{
  return (g->c[i] >> j) & 1u;
}



/* add an edge between i and j */
INLINE void dg_link(dg_t *g, int i, int j)
{
  g->c[i] |= (1u << j);
  g->c[j] |= (1u << i);
}



/* remove an edge between i and j */
INLINE void dg_unlink(dg_t *g, int i, int j)
{
  g->c[i] &= ~(1u << j);
  g->c[j] &= ~(1u << i);
}


/* construct `sg' by removing vertex `i0' from `g' */
INLINE dg_t *dg_shrink1(dg_t *sg, dg_t *g, int i0)
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
INLINE int dg_deg(const dg_t *g, int i)
{
  return bitcount(g->c[i]);
}


/* ascending */
int intcmp(const void *a, const void *b)
{
  return *(int *) b - *(int *) a;
}


/* get the degree sequence
 * return the number of partitions */
INLINE int dg_degseq(const dg_t *g, int *degseq)
{
  int i, sum = 0, n = g->n;

  for (i = 0; i < n; i++) {
    degseq[i] = bitcount(g->c[i]);
    sum += degseq[i];
  }
  qsort(degseq, n, sizeof(degseq[0]), intcmp);
  return sum / 2;
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
  int i;

  for (i = 0; i < g->n; i++)
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
  int i;
  code_t mask = (1 << g->n) - 1;

  for (i = 0; i < g->n; i++)
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
        dg_link(g, i, j);
      if (++ib == CODEBITS) ib = 0, c++;
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
  die_if (n >= CODEBITS, "do not support %d atoms\n", n);
  xnew(g->c, g->n);
  dg_empty(g);
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
  dg_t *a;

  a = dg_open(b->n);
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
INLINE int dg_connected(dg_t *g)
{
  code_t c = 1u, done = 1u, stack = 1u, mask = (1 << g->n) - 1;

  while (stack) { /* the next in the stack */
    int k = bitfirst(stack); /* first 1 bit */
    c |= g->c[k];
    done |= 1 << k; /* we have included k into the cluster */
    stack |= g->c[k]; /* add new connections */
    stack &= ~done; /* remove done indices */
    /* printf("k %d, stack 0x%x, done 0x%x\n", k, stack, done); */
    if (((stack | done) ^ mask) == 0) return 1;
  }
  return bitcount(c) == g->n;
}



/* check if the subdiagram without i is connected */
INLINE int dg_connectedi(const dg_t *g, int i)
{
  code_t c = 1u, done = 1u, stack = 1u, mask = ((code_t) 1u << g->n) - 1;

  if (i == 0) done = stack = c = 2u;
  mask &= ~(1 << i);
  while (stack) { /* the next in the stack */
    int k = bitfirst(stack); /* first 1 bit */
    c |= g->c[k];
    done |= 1 << k; /* we have included k into the cluster */
    stack |= g->c[k]; /* add new connections */
    stack &= (~done) & (~(1 << i)); /* remove done indices */
    /* printf("k %d, stack 0x%x, done 0x%x, mask 0x%x\n", k, stack, done, mask); */
    if (((stack | done) ^ mask) == 0) return 1;
  }
  c &= ~(1 << i); /* erase the ith bit */
  return bitcount(c) == g->n - 1;
}



/* check if diagram is biconnected, bitwise version */
INLINE int dg_biconnected(const dg_t *g)
{
  int i, n = g->n;
  static int bc4[64] = {
  /*   0-1 0-2 0-1,2 0-3 0-1,3 0-2,3  0-1,2,3 */
    0,  0,  0,  0,    0,  0,    0,     0, /* no edge among 1, 2, 3 */
    0,  0,  0,  0,    0,  0,    0,     0, /* only 1-2 */
    0,  0,  0,  0,    0,  0,    0,     0, /* only 1-3 */
    0,  0,  0,  0,    0,  0,    1,     1, /* 1-2, 1-3 */
    0,  0,  0,  0,    0,  0,    0,     0, /* only 2-3 */
    0,  0,  0,  0,    0,  1,    0,     1, /* 1-2, 2-3 */
    0,  0,  0,  1,    0,  0,    0,     1, /* 1-3, 2-3 */
    0,  0,  0,  1,    0,  1,    1,     1, /* 1-2, 1-3, 2-3 */
  };

  if (n < 2) return 1;
  else if (n == 2) return (g->c[0] >> 1) & 1u;
  else if (n == 3) return ((g->c[0] | (g->c[1] >> 2)) & 7u) == 7u;
  else if (n == 4) return bc4[ ((g->c[0] >> 1) | ((g->c[1] >> 2) << 3) | ((g->c[2] >> 3) << 5)) & 0x3f ];

  for (i = 0; i < n; i++)
    if ( !dg_connectedi(g, i) ) return 0;
  return 1;
}



typedef struct {
  int ng; /* number of unique diagrams */
  short *map; /* map a diagram to the unique diagram
                 `short' works for n <= 8 */
  code_t *first; /* index of the first unique diagram */
} dgmap_t;

/* static diagram map for n <= DGMAP_NMAX */
#define DGMAP_NMAX 8
dgmap_t dgmap_[DGMAP_NMAX];



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

  die_if (n > DGMAP_NMAX, "n %d is too large\n", n);
  if (m->ng > 0) return 0; /* already initialized */

  npr = n * (n - 1) / 2;
  ng = (code_t) 1u << npr;
  xnew(m->map, ng);
  for (c = 0; c < ng; c++) m->map[c] = -1;
  xnew(m->first, sz = 256);

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
    if (m->ng >= sz) xrenew(m->first, sz += 256);
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
  printf("n %d, %d unique diagrams\n", n, m->ng);
  return 0;
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



/* check if diagram is biconnected 
 * faster implementation */
INLINE int dg_biconnected_lookup(const dg_t *g)
{
  int k, n = g->n;
  code_t c;
  dgmap_t *m = dgmap_ + n;
  /* biconnectivity of unique diagrams */
  static int *bc[DGMAP_NMAX + 1] = {NULL};

  if (bc[n] == NULL) {
    dg_t *g1 = dg_open(n);

    dgmap_init(m, n); /* compute unique maps */
    xnew(bc[n], m->ng);
    for (k = 0; k < m->ng; k++) {
      c = m->first[k];
      dg_decode(g1, &c);
      bc[n][k] = dg_biconnected(g1);
    }
    dg_close(g1);
  }
  dg_encode(g, &c);
  k = m->map[c];
  return bc[n][k];
}


#endif
