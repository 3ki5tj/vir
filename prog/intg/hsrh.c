/* first few virial coefficients of the 3D hard spheres system
 * by direct Monte Carlo integration of Ree-Hoover diagrams
 * Copyright (C) 2013 Cheng Zhang */
/* #define `D' here to simulate dimensions other than D = 3 */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"



/* third-order virial coefficients */
double vir3(long nsteps)
{
  long t, acc = 0;
  rvn_t v, u;
  double x;

  for (t = 1; t <= nsteps; t++) {
    rvn_rndball0(u);
    rvn_inc(u, rvn_rndball0(v));
    if (rvn_sqr(u) < 1) acc++;
  }
  x = -2 * (-acc/6.) / nsteps;
  printf("vir3: %.8f, %.7f, %.9f, acc %g\n",
      x, x * pow(2, 2), x * pow(8, 2), 1. * acc / nsteps);
  return x;
}



/* fourth-order virial coefficients */
double vir4(long nsteps)
{
  long t, acc44 = 0, acc46 = 0;
  int i, lnk02, lnk13;
  rvn_t u, v[4];
  double x;

  for (t = 1; t <= nsteps; t++) {
    rvn_zero(v[0]);
    for (i = 0; i < 3; i++)
      rvn_add(v[i+1], v[i], rvn_rndball0(u));
    if (rvn_sqr(v[3]) < 1) {
      lnk02 = ( rvn_sqr(v[2]) < 1 );
      lnk13 = ( rvn_dist2(v[1], v[3]) < 1 );
      if (lnk02 && lnk13) acc46++;
      else if (!lnk02 && !lnk13) acc44++;
    }
  }
  printf("acc44 %g, acc46 %g\n", 1.*acc44/nsteps, 1.*acc46/nsteps);
  x = -3./24 * (3 * acc44 - 2 * acc46) / nsteps;
  printf("vir4: %.8f, %.7f, %.9f\n", x, x * pow(2, 3), x * pow(8, 3));
  return x;
}



/* fifth-order virial coefficients */
double vir5(long nsteps)
{
  rvn_t u, v[5];
  double x, y, tot, norm1, norm2;
  long t, acc[6] = {0};
  int i, mul[5] = {-12, 10, 60, 45, -6}, map[64];
  char *names[5] = {"A55", "A56", "A57", "A58", "A510"};
  unsigned code;

  for (i = 0; i < 64; i++)
    map[i] = 5; /* `5' means the un-biconnected diagrams */
  map[0x3f] = 4; /* fully connected */
  map[0x01] = 0; /* ring */
  map[0x14] = 1; /* house */

  /* one parallel connected lines in a ring; for each vertex,
   * the two lines form a cross from the opposite vertices */
  map[0x29] = 2;
  map[0x25] = 2;
  map[0x15] = 2;
  map[0x13] = 2;
  map[0x0b] = 2;

  /* two parallel wiggly lines */
  map[0x36] = 3;
  map[0x17] = 3;
  map[0x1b] = 3;
  map[0x2b] = 3;
  map[0x2d] = 3;
  map[0x35] = 3;

  rvn_zero(v[0]);
  for (t = 1; t <= nsteps; t++) {
    rvn_add(v[1], v[0], rvn_rndball0(u));
    rvn_add(v[2], v[1], rvn_rndball0(u));
    rvn_add(v[3], v[2], rvn_rndball0(u));
    rvn_add(v[4], v[3], rvn_rndball0(u));

#define CNT(i, j) (rvn_dist2(v[i], v[j]) < 1)
    code = CNT(0, 4)       | (CNT(0, 2) << 1) | (CNT(0, 3) << 2)
        | (CNT(1, 3) << 3) | (CNT(1, 4) << 4) | (CNT(2, 4) << 5);

    acc[ map[code] ]++;
  }
  acc[2] /= 5;
  acc[3] /= 6;

  norm1 = -4./120 * pow(2, 4);
  norm2 = -4./120 * pow(8, 4);
  for (tot = 0, i = 0; i < 5; i++) {
    x = 1. * acc[i] / nsteps;
    tot += y = x * mul[i];
    printf("%-5s acc%d: x %10.8f, mul %3d, x*mul %12.8f (%12.7f, "
        "%12.6f), tot %12.8f (%12.7f, %12.6f)\n",
        names[i], i, x, mul[i], y, norm1 * y, norm2 * y,
        tot, norm1 * tot, norm2 * tot);
  }
  tot *= -4./120;
  printf("vir5: %.8f, %.7f, %.6f\n", tot, tot * pow(2, 4), tot * pow(8, 4));
  return tot;
}



/* six-order virial coefficients */
double vir6(long nsteps)
{
#define NDG6 23 /* number of diagrams */
  long t, acc[NDG6 + 1] = {0};
  rvn_t u, v[6], vb[6];
  double x, y, tot, norm1, norm2;
  int i, j, k, l, ipr, npr, *map, *mapb;
  unsigned lnk, lnkb;
  struct {
    int sc, degen;
    int dup; /* occurrence in sampling */
    int npr; /* number of wiggly edges */
    int id[16][2];
    const char name[16];
    /* computed values */
    int mul;
    unsigned lnk; /* connectivity */
  } db[NDG6] = {
    /* cf. Reformulation of the Virial Series for Classical Fluids, Fig. 2
     * F. H. Ree and W. G. Hoover, J. Chem. Phys. 41, 6, 1635-1645
     * The diagrams have been reordered according to Table II. in
     * Virial Coefficients for hard spheres and hard disks,
     * F. H. Ree and W. G. Hoover, J. Chem. Phys. 40, 4 939-950 */
    {  24,   1,  1, 0, {}, "15"},
    { -12,  45,  1, 2, {{0, 4}, {1, 3}}, "13"},

    {  16,  15,  1, 3, {{0, 5}, {1, 3}, {2, 4}}, "12a"},
    {   8, 180,  1, 3, {{0, 2}, {2, 4}, {3, 5}}, "12b"},

    {  -4,  90,  1, 4, {{0, 2}, {2, 4}, {1, 3}, {3, 5}}, "11a"},
    {  -5, 180,  1, 4, {{0, 2}, {0, 4}, {1, 4}, {3, 5}}, "11b"},
    {  -4,  60,  1, 4, {{0, 2}, {2, 4}, {0, 4}, {3, 5}}, "11c"},
    {  -6,  60,  1, 4, {{0, 2}, {0, 3}, {0, 4}, {1, 5}}, "11d"},

    {  -1, 360,  1, 5, {{0, 3}, {3, 5}, {5, 2}, {2, 4}, {4, 1}}, "10a"},
    {   4,  45,  1, 5, {{0, 3}, {1, 4}, {4, 2}, {2, 5}, {5, 1}}, "10b"},
    {  -4,  72,  1, 5, {{0, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 0}}, "10c"}, /* 0-1-2-3-4-0 */
    {   3, 180,  1, 5, {{0, 2}, {0, 3}, {0, 5}, {3, 5}, {1, 4}}, "10d"},

    {   4,  10,  1, 6, {{0, 2}, {0, 4}, {2, 4}, {1, 3}, {3, 5}, {1, 5}}, "9a"},
    {   4,  60,  1, 6, {{0, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}, {5, 0}}, "9b"},
    {   1, 360,  1, 6, {{0, 3}, {0, 5}, {0, 2}, {2, 5}, {1, 3}, {1, 4}}, "9c"},
    {   3, 360,  1, 6, {{0, 5}, {0, 3}, {3, 1}, {1, 4}, {4, 2}, {2, 0}}, "9d"},
    {  -2,  90,  1, 6, {{0, 2}, {0, 4}, {0, 5}, {2, 4}, {2, 5}, {1, 3}}, "9e"},

    {  -2, 360,  1, 7, {{1, 3}, {3, 5}, {5, 2}, {2, 0}, {0, 4}, {4, 1}, {0, 5}}, "8a"},
    {  -1,  90,  1, 7, {{0, 2}, {2, 4}, {4, 5}, {1, 3}, {3, 5}, {5, 1}, {2, 5}}, "8b"},
    {  -2, 180,  1, 7, {{0, 4}, {0, 5}, {1, 3}, {1, 4}, {2, 4}, {2, 5}, {3, 5}}, "8c"},
    {   1,  15,  1, 7, {{0, 2}, {1, 3}, {1, 4}, {1, 5}, {3, 4}, {3, 5}, {4, 5}}, "8d"},
    /* the above diagram does not use the chain topology */

    {   1, 180,  1, 8, {{0, 2}, {0, 4}, {0, 5}, {1, 3}, {1, 4}, {2, 4}, {2, 5}, {3, 5}}, "7"},

    {   1,  60,  1, 9, {{0, 2}, {0, 3}, {0, 4}, {1, 3}, {1, 4}, {1, 5}, {2, 4}, {2, 5}, {3, 5}}, "6"},
  };

  /* alternative diagram structures with chain backbone */
  struct {
    int dgid; /* diagram id */
    int npr; /* number of edges */
    int id[16][2]; /* edges */
  } gchn[] = { /* negative edges */
    {1, 2, {{0, 5}, {1, 3}}},
    {1, 2, {{0, 5}, {1, 4}}},
    {1, 2, {{0, 5}, {2, 4}}},
    {1, 2, {{0, 2}, {1, 3}}},
    {1, 2, {{0, 2}, {1, 4}}},
    {1, 2, {{0, 2}, {1, 5}}},
    {1, 2, {{0, 2}, {3, 5}}},
    {1, 2, {{0, 3}, {1, 4}}},
    {1, 2, {{0, 3}, {1, 5}}},
    {1, 2, {{0, 3}, {2, 4}}},
    {1, 2, {{0, 3}, {2, 5}}},
    /* {1, 2, {{0, 4}, {1, 3}}}, */
    {1, 2, {{0, 4}, {1, 5}}},
    {1, 2, {{0, 4}, {2, 5}}},
    {1, 2, {{0, 4}, {3, 5}}},
    {1, 2, {{1, 3}, {2, 4}}},
    {1, 2, {{1, 3}, {2, 5}}},
    {1, 2, {{1, 4}, {2, 5}}},
    {1, 2, {{1, 4}, {3, 5}}},
    {1, 2, {{1, 5}, {2, 4}}},

    /* {2, 3, {{0, 5}, {1, 3}, {2, 4}}}, */
    {2, 3, {{0, 2}, {1, 4}, {3, 4}}},
    {2, 3, {{0, 3}, {1, 4}, {2, 5}}},
    {2, 3, {{0, 3}, {1, 5}, {2, 4}}},
    {2, 3, {{0, 4}, {2, 5}, {1, 3}}},

    /* one isolated wiggly line, and an isolated wiggly-line hat */
    {3, 3, {{0, 5}, {1, 3}, {1, 4}}},
    {3, 3, {{0, 5}, {1, 4}, {2, 4}}},
    {3, 3, {{0, 2}, {1, 3}, {1, 4}}},
    {3, 3, {{0, 2}, {1, 3}, {1, 5}}},
    {3, 3, {{0, 2}, {1, 4}, {1, 5}}},
    {3, 3, {{0, 2}, {1, 3}, {3, 5}}},
    {3, 3, {{0, 2}, {1, 5}, {3, 5}}},
    {3, 3, {{0, 3}, {1, 4}, {1, 5}}},
    {3, 3, {{0, 3}, {2, 4}, {2, 5}}},
    {3, 3, {{0, 3}, {4, 1}, {4, 2}}},
    {3, 3, {{0, 3}, {5, 1}, {5, 2}}},
    {3, 3, {{0, 4}, {1, 3}, {1, 5}}},
    {3, 3, {{0, 4}, {1, 3}, {3, 5}}},
    {3, 3, {{0, 4}, {1, 5}, {3, 5}}},
    {3, 3, {{1, 3}, {0, 2}, {0, 4}}},
    {3, 3, {{1, 3}, {0, 2}, {2, 4}}},
    {3, 3, {{1, 3}, {0, 2}, {2, 5}}},
    {3, 3, {{1, 3}, {2, 4}, {2, 5}}},
    {3, 3, {{1, 3}, {0, 4}, {2, 4}}},
    {3, 3, {{1, 3}, {0, 2}, {0, 5}}},
    {3, 3, {{1, 3}, {0, 4}, {0, 5}}},
    {3, 3, {{1, 3}, {0, 5}, {2, 5}}},
    {3, 3, {{1, 4}, {0, 2}, {0, 3}}},
    {3, 3, {{1, 4}, {0, 2}, {2, 5}}},
    {3, 3, {{1, 4}, {0, 3}, {3, 5}}},
    {3, 3, {{1, 4}, {2, 5}, {3, 5}}},
    {3, 3, {{1, 4}, {0, 2}, {0, 5}}},
    {3, 3, {{1, 4}, {0, 3}, {0, 5}}},
    {3, 3, {{1, 4}, {0, 5}, {2, 5}}},
    {3, 3, {{1, 4}, {0, 5}, {3, 5}}},
    {3, 3, {{1, 5}, {0, 2}, {0, 3}}},
    {3, 3, {{1, 5}, {0, 3}, {0, 4}}},
    {3, 3, {{1, 5}, {0, 2}, {0, 4}}},
    {3, 3, {{1, 5}, {0, 2}, {2, 4}}},
    {3, 3, {{1, 5}, {0, 4}, {2, 4}}},
    {3, 3, {{2, 4}, {1, 3}, {1, 5}}},
    {3, 3, {{2, 4}, {1, 3}, {3, 5}}},
    {3, 3, {{2, 4}, {0, 3}, {3, 5}}},
    {3, 3, {{2, 4}, {0, 3}, {1, 3}}},
    {3, 3, {{2, 4}, {1, 5}, {3, 5}}},
    {3, 3, {{2, 4}, {0, 5}, {0, 3}}},
    {3, 3, {{2, 4}, {0, 5}, {1, 5}}},
    {3, 3, {{2, 4}, {0, 5}, {3, 5}}},
    {3, 3, {{2, 5}, {0, 3}, {0, 4}}},
    {3, 3, {{2, 5}, {1, 3}, {1, 4}}},
    {3, 3, {{2, 5}, {0, 3}, {1, 3}}},
    {3, 3, {{2, 5}, {0, 4}, {1, 4}}},
    {3, 3, {{3, 5}, {0, 2}, {0, 4}}},
    /* {3, 3, {{3, 5}, {0, 2}, {2, 4}}}, */
    {3, 3, {{3, 5}, {0, 4}, {1, 4}}},
    {3, 3, {{3, 5}, {0, 4}, {2, 4}}},
    {3, 3, {{3, 5}, {1, 4}, {2, 4}}},

    /* {7, 4, {{0, 2}, {0, 3}, {0, 4}, {1, 5}}}, */
    {7, 4, {{0, 2}, {0, 4}, {0, 5}, {1, 3}}},
    {7, 4, {{0, 2}, {0, 3}, {0, 5}, {1, 4}}},
    {7, 4, {{5, 3}, {5, 2}, {5, 1}, {4, 0}}},
    {7, 4, {{5, 3}, {5, 1}, {5, 0}, {4, 2}}},
    {7, 4, {{5, 3}, {5, 2}, {5, 0}, {4, 1}}},
    {7, 4, {{1, 3}, {1, 4}, {1, 5}, {0, 2}}},
    {7, 4, {{4, 2}, {4, 1}, {4, 0}, {5, 3}}},
    {7, 4, {{2, 0}, {2, 4}, {2, 5}, {1, 3}}},
    {7, 4, {{3, 5}, {3, 1}, {3, 0}, {4, 2}}},

    {10, 5, {{0, 3}, {3, 5}, {5, 2}, {2, 4}, {4, 0}}}, /* 0-2-3-4-5-0 */
    {10, 5, {{0, 3}, {3, 5}, {5, 1}, {1, 4}, {4, 0}}}, /* 0-1-3-4-5-0 */
    {10, 5, {{0, 2}, {2, 5}, {5, 1}, {1, 4}, {4, 0}}}, /* 0-1-2-4-5-0 */
    {10, 5, {{0, 2}, {2, 5}, {5, 1}, {1, 3}, {3, 0}}}, /* 0-1-2-3-5-0 */
    {10, 5, {{5, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}}}, /* 5-1-2-3-4-5 */
    {10, 5, {{0, 5}, {5, 1}, {1, 4}, {4, 2}, {2, 0}}}, /* 0-5-1-4-2-0 */
    {10, 5, {{5, 0}, {0, 4}, {4, 1}, {1, 3}, {3, 5}}}, /* 5-0-4-1-3-5 */

    {15, 6, {{0, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 0}, {1, 5}}}, /* 0-1-2-3-4-0, 5 */
    {15, 6, {{0, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 0}, {2, 5}}}, /* 0-1-2-3-4-0, 5 */
    {15, 6, {{0, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 0}, {3, 5}}}, /* 0-1-2-3-4-0, 5 */
    {15, 6, {{0, 3}, {3, 5}, {5, 2}, {2, 4}, {4, 0}, {3, 1}}}, /* 0-2-3-4-5-0, 1 */
    {15, 6, {{0, 3}, {3, 5}, {5, 2}, {2, 4}, {4, 0}, {4, 1}}}, /* 0-2-3-4-5-0, 1 */
    {15, 6, {{0, 3}, {3, 5}, {5, 2}, {2, 4}, {4, 0}, {5, 1}}}, /* 0-2-3-4-5-0, 1 */
    {15, 6, {{0, 3}, {3, 5}, {5, 1}, {1, 4}, {4, 0}, {0, 2}}}, /* 0-1-3-4-5-0, 2 */
    {15, 6, {{0, 3}, {3, 5}, {5, 1}, {1, 4}, {4, 0}, {4, 2}}}, /* 0-1-3-4-5-0, 2 */
    {15, 6, {{0, 3}, {3, 5}, {5, 1}, {1, 4}, {4, 0}, {5, 2}}}, /* 0-1-3-4-5-0, 2 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 4}, {4, 0}, {0, 3}}}, /* 0-1-2-4-5-0, 3 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 4}, {4, 0}, {1, 3}}}, /* 0-1-2-4-5-0, 3 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 4}, {4, 0}, {5, 3}}}, /* 0-1-2-4-5-0, 3 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 3}, {3, 0}, {0, 4}}}, /* 0-1-2-3-5-0, 4 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 3}, {3, 0}, {1, 4}}}, /* 0-1-2-3-5-0, 4 */
    {15, 6, {{0, 2}, {2, 5}, {5, 1}, {1, 3}, {3, 0}, {2, 4}}}, /* 0-1-2-3-5-0, 4 */
    {15, 6, {{5, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}, {2, 0}}}, /* 5-1-2-3-4-5, 0 */
    {15, 6, {{5, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}, {3, 0}}}, /* 5-1-2-3-4-5, 0 */
    {15, 6, {{5, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}, {4, 0}}}, /* 5-1-2-3-4-5, 0 */
    {15, 6, {{5, 2}, {2, 4}, {4, 1}, {1, 3}, {3, 5}, {5, 0}}}, /* 5-1-2-3-4-5, 0 */
    {15, 6, {{0, 5}, {5, 1}, {1, 4}, {4, 2}, {2, 0}, {3, 0}}}, /* 0-5-1-4-2-0, 3 */
    {15, 6, {{0, 5}, {5, 1}, {1, 4}, {4, 2}, {2, 0}, {3, 1}}}, /* 0-5-1-4-2-0, 3 */
    {15, 6, {{0, 5}, {5, 1}, {1, 4}, {4, 2}, {2, 0}, {3, 5}}}, /* 0-5-1-4-2-0, 3 */
    {15, 6, {{5, 0}, {0, 4}, {4, 1}, {1, 3}, {3, 5}, {2, 0}}}, /* 5-0-4-1-3-5, 2 */
    {15, 6, {{5, 0}, {0, 4}, {4, 1}, {1, 3}, {3, 5}, {2, 4}}}, /* 5-0-4-1-3-5, 2 */
    {15, 6, {{5, 0}, {0, 4}, {4, 1}, {1, 3}, {3, 5}, {2, 5}}}, /* 5-0-4-1-3-5, 2 */

    {-1, 0, {}},
  }, gchp[] = {
    {19, 3, {{0, 5}, {0, 2}, {1, 3}}},
    {19, 3, {{0, 5}, {1, 3}, {2, 4}}},
    {19, 3, {{0, 5}, {2, 4}, {3, 5}}},
    {19, 3, {{0, 5}, {3, 5}, {0, 4}}},
    {19, 3, {{0, 5}, {0, 4}, {1, 5}}},
    {21, 2, {{2, 5}, {0, 4}}},
    {-1, 0, {}},
  };

  npr = (1 << 15);
  xnew(map, npr);
  xnew(mapb, npr);
  for (ipr = 0; ipr < npr; ipr++)
    map[ipr] = mapb[ipr] = NDG6; /* invalid entry */

  /* compute multiplicity */
  for (k = 0; k < NDG6; k++) {
    db[k].mul = db[k].sc * db[k].degen * ((15 - db[k].npr) % 2 ? -1 : 1);
    db[k].lnk = (1u << 15) - 1; /* fully connected diagram */
    for (ipr = 0; ipr < db[k].npr; ipr++) {
      i = getpairindex( db[k].id[ipr][0], db[k].id[ipr][1], 6 );
      db[k].lnk &= ~(1u << i); /* unlink the bond */
    }
    if (k == 20) mapb[ db[k].lnk ] = k; /* make the code to the graph */
    else map[ db[k].lnk ] = k;
  }

  /* add alternative diagrams (specified by wiggly lines) */
  for (l = 0; (k = gchn[l].dgid) >= 0; l++) {
    lnk = 0x7fff;
    for (ipr = 0; ipr < gchn[l].npr; ipr++) {
      i = getpairindex( gchn[l].id[ipr][0], gchn[l].id[ipr][1], 6 );
      lnk &= ~(1u << i); /* unlink the bond */
    }
    map[ lnk ] = k;
    db[k].dup++;
  }

  /* add alternative diagrams (specified by wiggly lines) */
  for (l = 0; (k = gchp[l].dgid) >= 0; l++) {
    lnk = 0;
    for (j = 0; j < 5; j++) { /* construct the chain */
      i = getpairindex(j, j + 1, 6 );
      lnk |= 1u << i; /* link the bond */
    }
    for (ipr = 0; ipr < gchp[l].npr; ipr++) {
      i = getpairindex( gchp[l].id[ipr][0], gchp[l].id[ipr][1], 6 );
      lnk |= 1u << i; /* unlink the bond */
    }
    map[ lnk ] = k;
    db[k].dup++;
  }

  for (k = 0; k < NDG6; k++)
    printf("%2d A6_%-4s lnk 0x%04x, mul %+5d, dup %3d\n",
        k, db[k].name, db[k].lnk, db[k].mul, db[k].dup);

  rvn_zero(v[0]);
  rvn_zero(vb[0]);
  for (t = 1; t <= nsteps; t++) {
    /* spanning tree of chain topology */
    for (i = 0; i < 5; i++)
      rvn_add(v[i + 1], v[i], rvn_rndball0(u));
    rvn_copy(vb[1], v[1]);
    rvn_copy(vb[2], v[2]);
    rvn_copy(vb[3], v[3]);
    rvn_add(vb[4], vb[2], rvn_rndball0(u));
    rvn_add(vb[5], vb[2], rvn_rndball0(u));

    /* encode the connectivity */
    for (lnk = lnkb = 0, ipr = 0, i = 0; i < 5; i++)
      for (j = i + 1; j < 6; j++) {
        if (rvn_dist2(v[i], v[j]) < 1)
          lnk |= 1 << ipr;
        if (rvn_dist2(vb[i], vb[j]) < 1)
          lnkb |= 1 << ipr;
        ipr++;
      }

    acc[ map[ lnk ] ]++;
    acc[ mapb[ lnkb ] ]++;
#if 0
    /* fully-connected diagram */
    if (lnkb == db[0].lnk) db[0].acc++;
    if (lnk == db[0].lnk) db[0].acc++;

    for (k = 1; k < NDG6; k++) /* accumulate acceptance rates */
      if (k == 20) {
        if (db[k].lnk == lnkb)
          db[k].acc++;
      } else if (db[k].lnk == lnk) {
        db[k].acc++;
      }
#endif
  }

  norm1 = -5./720 * pow(2, 5);
  norm2 = -5./720 * pow(8, 5);
  for (tot = 0, k = 0; k < NDG6; k++) {
    x = 1. * acc[k] / nsteps / db[k].dup;
    tot += y = x * db[k].mul;
    printf("%2d A6_%-3s: acc %10.8f, sc %+3d, sc*dg %+5d mul %+5d "
        "x*mul %12.8f (%12.7f, %12.6f), tot %12.8f (%12.7f, %12.6f)\n",
        k, db[k].name, x, db[k].sc, db[k].sc * db[k].degen,
        db[k].mul, y, norm1 * y, norm2 * y,
        tot, norm1 * tot, norm2 * tot);
  }
  tot *= -5./720;
  printf("vir6: %.8f, %.7f, %.6f\n", tot, tot * pow(2, 5), tot * pow(8, 5));
  free(map); free(mapb);
  return tot;
}



int main(int argc, char *argv[])
{
  long nsteps = 10000000;
  int n = 6;

  /* Metropolis, Rosenbluth, Roksenbluth, Teller & Teller, J. Chem. Phys. 21, 6, 1087-1092, 1953
   * Ree, Hoover, J. Chem. Phys. 46, 4, 939-950, 1964
   * http://www.sklogwiki.org/SklogWiki/index.php/Hard_sphere:_virial_coefficients
   * B2  = b0 = 2 pi/3
   * B3  = 5/8 b0^2
   * B4  = 0.2869495 b0^3
   * B5  = 0.110252(1) b0^4
   * B6  = 0.03888198(91) b0^5
   * B7  = 0.01302354(91) b0^6
   * B8  = 0.0041832(11) b0^7
   * B9  = 0.0013094(13) b0^8
   * B10 = 0.0004035(15) b0^9 */
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nsteps = strtol(argv[2], NULL, 10);
  printf("n = %d, %ld Monte Carlo moves\n", n, nsteps);

  if (n == 3) vir3(nsteps);
  else if (n == 4) vir4(nsteps);
  else if (n == 5) vir5(nsteps);
  else if (n == 6) vir6(nsteps);
  mtsave(NULL);
  return 0;
}

