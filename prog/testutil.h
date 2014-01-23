#ifndef TESTUTIL_H__
#define TESTUTIL_H__
/* utilities for test programs */
#include "dg.h"



#define dg_rndswitchedge(g, ned, nedmax, rnp) \
  dg_rndswitchedge0(g, ned, 0, 1, nedmax, rnp)

INLINE int dg_rndswitchedge0(dg_t *g, int *ned,
    int nedmin, double rnm, int nedmax, double rnp)
{
  int i = 0, j = 0, n = g->n, ipr, npr, acc;

  npr = n * (n - 1) / 2;
  /* randomly switch an edge */
  ipr = (int) (npr * rnd0());
  parsepairindex(ipr, n, &i, &j);
  if (dg_linked(g, i, j)) { /* unlink (i, j) */
    dg_unlink(g, i, j);
    acc = dg_biconnected(g);
    if ( acc && *ned <= nedmin )
      acc = (rnd0() < rnm);
    if ( acc ) {
      (*ned)--;
    } else { /* link back the edge */
      dg_link(g, i, j);
    }
  } else { /* link (i, j) */
    acc = 1;
    if (*ned >= nedmax) /* avoid increasing edges */
      acc = (rnd0() < rnp);
    if (acc) {
      (*ned)++;
      dg_link(g, i, j);
    }
  }
  return acc;
}



INLINE void dg_linkpairs(dg_t *g, int pair[][2])
{
  int i;

  for (i = 0; pair[i][0] >= 0; i++)
    dg_link(g, pair[i][0], pair[i][1]);
}



/* special ring content formulas
 *
 * a~b*k
 * = (n - 1)!/2*1F1(-k, -n+1, -2)
 *
 * a~b~c      d~e*k
 * = (n - 1)!/2*1F1(-k - 2, -n + 1; -2)
 * - (n - 3)!*1F1(-k, -n + 3; -2)
 *
 * a~b~c~d    e~f*k
 * = (n - 1)!/2*1F1(-k - 3, -n + 1; -2)
 * - 2*(n - 3)!*1F1(-k, -n + 3; -2)
 * + 3*(n - 4)!*1F1(-k, -n + 4; -2)
 *
 * a~b~c~d    e~f*k
 * = (n - 1)!/2*1F1(-k - 3, -n + 1; -2)
 * - 2*(n - 3)!*1F1(-k, -n + 3; -2)
 * + 3*(n - 4)!*1F1(-k, -n + 4; -2)
 *
 * a~b~c~d~e  f~g*k
 * =1/2*(n - 1)!*1F1(-k - 4, -n + 1; -2)
 * -  3*(n - 3)!*1F1(-k, -n + 3; -2)
 * + 10*(n - 4)!*1F1(-k, -n + 4, -2)
 * -  7*(n - 5)!*1F1(-k, -n + 5, -2)
 *
 * a~b~c~d~e~f g~h*k
 * =1/2*(n - 1)!*1F1(-k-5, -n + 1; -2)
 * -  4*(n - 3)!*1F1(-k, -n + 3, -2)
 * + 21*(n - 4)!*1F1(-k, -n + 4, -2)
 * - 32*(n - 5)!*1F1(-k, -n + 5, -2)
 * + 15*(n - 6)!*1F1(-k, -n + 6, -2)
 *
 * k-clique
 *  1/2 (n - 1 - k)! (n - k)! / (n - 2k)!
 *
 * */


#define DGREF_NPRMAX 64

typedef struct {
  int npr; /* number of wiggly lines */
  int id[DGREF_NPRMAX][2]; /* vertex pairs in wiggly lines */
  int fb; /* the correct signed star content */
  int nr; /* the correct ring content */
  int cs; /* contains a clique separator  */
} dgref_t;

/* test sets */
dgref_t dgref4[] = {
  {0, {{0, 0}}, -2, 3},
  {-1, {{0, 1}}, 0, 1, 1},
  {-2, {{0, 1}, {2, 3}}, 1, 1},
  {-2, {{1, 3}, {2, 0}}, 1, 1},
  {DGREF_NPRMAX, {{0, 0}}, 0, 0, 0},
};

dgref_t dgref5[] = {
  {0, {{0, 0}}, -6, 12},
  {-1, {{0, 1}}, 0, 6, 1},
  {-2, {{0, 1}, {2, 3}}, 3, 4}, /* a~b c~d */
  {-2, {{0, 1}, {1, 3}}, 0, 2, 1}, /* a~b~c */
  {-2, {{0, 4}, {2, 4}}, 0, 2, 1}, /* a~b~c */
  {-3, {{0, 1}, {2, 3}, {3, 4}}, 2, 2}, /* a~b~c d~e */
  {-3, {{0, 1}, {1, 3}, {2, 4}}, 2, 2}, /* a~b~c d~e */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, 0, 1, 1}, /* 4-ring */
  {6, {{0, 2}, {2, 3}, {3, 1}, {1, 4}, {4, 0}, {0, 3}}, 0, 1, 1},
  {6, {{0, 1}, {1, 3}, {3, 4}, {4, 0}, {1, 2}, {2, 4}}, 1, 0},
  {5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -1, 1},
  {DGREF_NPRMAX, {{0, 0}}, 0, 0},
};

dgref_t dgref6[] = {
  {0, {{0, 0}}, -24, 60},
  {-1, {{0, 1}}, 0, 36, 1}, /* a~b */
  {-2, {{0, 1}, {2, 3}}, 12, 24}, /* a~b c~d */
  {-2, {{0, 4}, {2, 4}}, 0, 18, 1}, /* a~b~c */
  {-3, {{0, 1}, {2, 3}, {4, 5}}, 16, 16}, /* a~b c~d e~f */
  {-3, {{0, 1}, {2, 3}, {3, 4}}, 8, 14}, /* a~b~c d~e */
  {-3, {{0, 1}, {1, 2}, {0, 2}}, 0, 6, 1}, /* 3-clique */
  {-4, {{0, 1}, {1, 2}, {0, 2}, {3, 4}}, 4, 6}, /* 3-clique a~b */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {4, 5}}, 5, 8}, /* a~b~c~d e~f */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {3, 1}}, 0, 2, 1}, /* 4-ring */
  {-5, {{0, 1}, {1, 2}, {0, 2}, {3, 4}, {3, 5}}, 0, 6, 0}, /* 3-clique a~b~c */
  {-5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, -4, 5}, /* 5-ring */
  {-6, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 0, 0, 1}, /* 4-clique */
  {-6, {{0, 1}, {0, 2}, {1, 2}, {3, 4}, {4, 5}, {5, 3}}, -4, 6}, /* two 3-rings */
  {-7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {2, 5}, {3, 5}}, 0, 0, 1},
  {7, {{4, 2}, {2, 3}, {3, 5}, {5, 1}, {1, 0}, {0, 4}, {1, 3}}, 0, 1, 1},
  {7, {{4, 3}, {3, 5}, {5, 1}, {1, 2}, {2, 4}, {0, 3}, {0, 1}}, -1, 0},
  {6, {{4, 5}, {5, 1}, {1, 3}, {3, 2}, {2, 0}, {0, 4}}, 1, 1},
  {DGREF_NPRMAX, {{0, 0}}, 0, 0},
};

dgref_t dgref7[] = {
  {0, {{0, 0}}, -120, 360},
  {-1, {{0, 1}}, 0, 240, 1}, /* a~b */
  {-2, {{0, 1}, {2, 3}}, 60, 168}, /* a~b c~d */
  {-2, {{4, 5}, {5, 6}}, 0, 144, 1}, /* a~b~c */
  {-3, {{0, 1}, {1, 2}, {3, 4}}, 40, 108}, /* a~b~c d~e */
  {-3, {{0, 1}, {2, 3}, {4, 5}}, 80, 120}, /* a~b c~d e~f */
  {-3, {{0, 1}, {1, 2}, {0, 2}}, 0, 72, 1}, /* 3-clique of e-bonds */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {4, 5}}, 25, 70}, /* a~b~c~d e~f */
  {-5, {{6, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 5}}, 14, 48}, /* a~b~c~d~e f~g */
  {-6, {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}, 0, 0, 1}, /* 3-clique */
  {-6, {{0, 1}, {0, 2}, {1, 2}, {3, 4}, {4, 5}, {5, 3}}, -20, 36}, /* two 3-rings */
  {-7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {2, 5}, {3, 5}}, 0, 8, 1},
  {-8, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {2, 5}, {3, 5}, {2, 4}}, 0, 0, 1},
  {8, {{4, 5}, {5, 3}, {3, 1}, {1, 2}, {2, 6}, {2, 0}, {0, 4}, {4, 6}}, 1, 0},
  {8, {{4, 5}, {5, 3}, {3, 1}, {1, 2}, {2, 6}, {6, 0}, {0, 4}, {4, 6}}, 0, 1, 1},
  {8, {{4, 5}, {5, 3}, {3, 1}, {1, 2}, {2, 6}, {6, 0}, {0, 4}, {2, 0}}, 0, 1, 1},
  {7, {{4, 2}, {2, 3}, {3, 5}, {5, 1}, {1, 6}, {6, 0}, {0, 4}}, -1, 1},
  {DGREF_NPRMAX, {{0, 0}}, 0, 0},
};

dgref_t dgref8[] = {
  {0, {{0, 0}}, -720, 2520},
  {-1, {{0, 1}}, 0, 1800, 1}, /* a~b */
  {-2, {{5, 6}, {4, 7}}, 360, 1320}, /* a~b b~c */
  {-2, {{4, 5}, {5, 6}}, 0, 1200, 1}, /* a~b~c */
  {-3, {{0, 1}, {1, 2}, {4, 5}}, 240, 912}, /* a~b~c d~e */
  {-3, {{0, 1}, {2, 3}, {4, 5}}, 480, 984}, /* a~b c~d e~f */
  {-4, {{0, 1}, {2, 3}, {4, 5}, {6, 7}}, 450, 744}, /* a~b c~d e~f g~h */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {4, 5}}, 150, 636}, /* a~b~c~d e~f */
  {-5, {{6, 1}, {1, 2}, {2, 3}, {4, 5}, {0, 7}}, 168, 496}, /* a~b~c~d e~f g~h */
  {-5, {{6, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 7}}, 84, 458}, /* a~b~c~d~e f~g */
  {-6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {6, 7}}, 13, 340}, /* a~b~c~d~e~f g~h */
  {9, {{5, 2}, {2, 4}, {4, 3}, {3, 1}, {1, 6}, {6, 0}, {0, 7}, {0, 5}, {2, 7}}, -1, 0},
  {9, {{5, 2}, {2, 4}, {4, 3}, {3, 1}, {1, 6}, {6, 0}, {0, 7}, {7, 5}, {1, 7}}, 0, 1, 1},
  {9, {{5, 2}, {2, 4}, {4, 3}, {3, 1}, {1, 6}, {6, 0}, {0, 7}, {7, 5}, {3, 6}}, 0, 1, 1},
  {8, {{7, 2}, {2, 3}, {3, 5}, {5, 1}, {1, 6}, {6, 0}, {0, 4}, {4, 7}}, 1, 1},
  {DGREF_NPRMAX, {{0, 0}}, 0, 0},
};

dgref_t dgref9[] = {
  {0, {{0, 0}}, -5040, 20160},
  {-1, {{3, 8}}, 0, 15120, 1}, /* a~b */
  {-2, {{5, 6}, {4, 7}}, 2520, 11520}, /* a~b c~d */
  {-2, {{4, 5}, {5, 6}}, 0, 10800, 1}, /* a~b~c */
  {-3, {{0, 1}, {1, 2}, {4, 5}}, 1680, 8400}, /* a~b~c d~e */
  {-3, {{0, 1}, {2, 3}, {4, 5}}, 3360, 8880}, /* a~b c~d e~f */
  {-4, {{4, 5}, {5, 6}, {1, 2}, {7, 8}}, 2100, 6576}, /* a~b~c d~e f~g */
  {-4, {{0, 1}, {2, 3}, {4, 5}, {6, 7}}, 3150, 6912}, /* a~b c~d e~f g~h */
  {-4, {{0, 1}, {1, 2}, {2, 3}, {4, 5}}, 1050, 6168}, /* a~b~c~d e~f */
  {-5, {{6, 1}, {1, 2}, {2, 3}, {4, 5}, {0, 7}}, 1176, 4896}, /* a~b~c~d e~f g~h */
  {-5, {{6, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 7}}, 588, 4620}, /* a~b~c~d~e f~g */
  {-6, {{6, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 7}, {5, 8}}, 518, 3704}, /* a~b~c~d~e f~g h~i */
  {-6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {6, 7}}, 91, 3526}, /* a~b~c~d~e~f g~h */
  {DGREF_NPRMAX, {{0, 0}}, 0, 0},
};

#define DGREF_NMIN 4
#define DGREF_NMAX 9

dgref_t *dgrefs[] = {NULL, NULL, NULL, NULL,
  dgref4, dgref5, dgref6, dgref7, dgref8, dgref9};



/* build a graph according to dgref_t */
INLINE void dgref_build(dg_t *g, dgref_t *ref)
{
  int j, npr = ref->npr;

  if (npr == DGREF_NPRMAX) return;
  /* build the graph */
  if ( npr <= 0 ) { /* e-bonds, removed edges */
    npr = -npr;
    dg_full(g);
    for (j = 0; j < npr; j++)
      dg_unlink(g, ref->id[j][0], ref->id[j][1]);
  } else { /* f-bonds needed edges */
    dg_empty(g);
    for (j = 0; j < npr; j++)
      dg_link(g, ref->id[j][0], ref->id[j][1]);
  }
}



/* adjust the rate of adding edges */
INLINE void adjustrnp(int ned, int nedmax, int t, int nstadj,
    int *good, int *tot, double *rnp)
{
  *good += (ned <= nedmax);
  *tot += 1;
  if (t % nstadj == 0 && 1. * (*good) / (*tot) < 0.5) { /* adjust rnp */
    *rnp *= 0.5;
    printf("adjusting rnp to %g, ned %d, good %g%%\n",
        *rnp, ned, 100.*(*good)/(*tot));
    *good = *tot = 0;
  }
}


/* adjust the rate of removing edges */
INLINE void adjustrnm(int ned, int nedmin, int t, int nstadj,
    int *good, int *tot, double *rnm)
{
  *good += (ned >= nedmin);
  *tot += 1;
  if (t % nstadj == 0 && 1. * (*good) / (*tot) < 0.5) { /* adjust rnp */
    *rnm *= 0.5;
    printf("adjusting rnm to %g, ned %d, good %g%%\n",
        *rnm, ned, 100.*(*good)/(*tot));
    *good = *tot = 0;
  }
}
#endif
