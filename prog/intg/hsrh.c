/* first few virial coefficients of the 3D hard spheres system
 * by direct Monte Carlo integration of Ree-Hoover diagrams
 * Copyright (C) 2013 Cheng Zhang */



/* #define `D' here to simulate dimensions other than D = 3 */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"



/* third-order virial coefficients */
static double vir3(long nsteps)
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
static double vir4(long nsteps)
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
static double vir5(long nsteps)
{
  rvn_t u, v[5];
  double x, y, tot, norm1, norm2;
  long t, acc[6] = {0}; /* acc[5] is reserved for invalid diagrams */
  int sc[5] = {-6, 3, -2, 1, 1}, degen[5] = {1, 15, -30, 10, -12};
  int i, mul[5], dup[5] = {1, 6, 7, 1, 1}, map[64];
  const char *names[5] = {"A5_10", "A5_8", "A5_7", "A5_6", "A5_5"};
  unsigned code;

  for (i = 0; i < 5; i++)
    mul[i] = sc[i] * degen[i];

  for (i = 0; i < 64; i++)
    map[i] = 5; /* `5' means invalid or un-biconnected diagrams */
  map[0x3f] = 0; /* fully connected */
  map[0x01] = 4; /* ring 0-4 */
  map[0x14] = 3; /* house 0-3, 1-4 */

  /* one wiggly line (L), one wiggly angle (A) */
  /* edge: ij (hex): 04(1) 02(2) 03(4) 13(8) 14(10) 24(20) */
  map[0x25] = 2; /*     #     L     #     A     A      #   */
  map[0x0b] = 2; /*     #     #     L     #     A      A   */
  map[0x15] = 2; /*     #     A     #     L     #      A   */
  map[0x34] = 2; /*     A     A     #     L     #      #   */
  map[0x16] = 2; /*     A     #     #     L     #      A   */
  map[0x29] = 2; /*     #     A     A     #     L      #   */
  map[0x13] = 2; /*     #     #     A     A     #      L   */

  /* two parallel wiggly lines */
  /* edge: ij (hex): 04(1) 02(2) 03(4) 13(8) 14(10) 24(20) */
  map[0x36] = 1; /*           #     #           #      #   */
  map[0x35] = 1; /*     #           #           #      #   */
  map[0x2d] = 1; /*     #           #     #            #   */
  map[0x2b] = 1; /*     #     #           #            #   */
  map[0x1b] = 1; /*     #     #           #     #          */
  map[0x17] = 1; /*     #     #     #           #          */

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

  norm1 = -4./120 * pow(2, 4);
  norm2 = -4./120 * pow(8, 4);
  for (tot = 0, i = 0; i < 5; i++) {
    x = 1. * acc[i] / nsteps / dup[i];
    tot += y = x * mul[i];

    printf("%2d %-5s acc %9.7f dup %d sc %+2d mul %+3d "
        "x*mul %10.7f (%9.6f,%10.5f) tot %10.7f (%9.6f,%10.5f)\n",
        i, names[i], x, dup[i], sc[i], mul[i], y, norm1 * y, norm2 * y,
        tot, norm1 * tot, norm2 * tot);
  }
  tot *= -4./120;
  printf("vir5: %.8f, %.7f, %.6f\n", tot, tot * pow(2, 4), tot * pow(8, 4));
  return tot;
}



#define NMAX 8
#define NPRMAX (NMAX*(NMAX-1)/2)

static unsigned encode(int n, int npr, int id[][2])
{
  int ipr;
  unsigned code = 0;

  for (ipr = 0; ipr < npr; ipr++)
    code |= 1 << getpairindex(id[ipr][0], id[ipr][1], n);
  return code;
}



/* given a backbone structure `bb', find the number of permutations
 * of the diagram of wiggly lines `pr' compatible with it */
static int mktop(int n, int npr, int pr[][2], int k, int map[],
    int nbb, int bb[][2], int verbose)
{
  int st[NMAX], used[NMAX] = {0}, top = 0, i, ipr, idn[NPRMAX][2], dup = 0;
  unsigned mask = 0, code;

  mask = encode(n, nbb, bb);
  for (i = 0; i < n; i++) st[i] = -1;

  /* iterate over all permutations of the diagram `pr'
   * see if it is compatible with the backbone `bb' */
  while (1) {
    if (top >= n) {
      for (ipr = 0; ipr < npr; ipr++) {
        idn[ipr][0] = st[ pr[ipr][0] ];
        idn[ipr][1] = st[ pr[ipr][1] ];
      }
      code = encode(n, npr, idn);
      if (map[code] >= 0 && map[code] != k) {
        printf("conflict for k %d, code %#x\n", k, code);
        for (ipr = 0; ipr < npr; ipr++)
          printf("(%d, %d) - (%d, %d)\n", pr[ipr][0], pr[ipr][1], idn[ipr][0], idn[ipr][1]);
        exit(1);
      }
      /* `mask' includes needed edges, `code' includes excluded edges */
      if ((mask & code) == 0 && map[code] < 0) {
        map[code] = k;
        dup++;
        if (verbose) {
          printf("k %2d, dup %2d: ", k, dup);
          for (ipr = 0; ipr < npr; ipr++)
            printf("{%d, %d}, ", idn[ipr][0], idn[ipr][1]);
          printf("\n");
        }
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
  if (verbose) printf("\n");
  return dup;
}



/* six-order virial coefficients */
static double vir6(long nsteps)
{
#define NDG6 23 /* number of diagrams */
  long t, acc[NDG6 + 1] = {0};
  rvn_t u, v[6], vb[6];
  double x, y, tot, norm1, norm2;
  int i, j, k, ipr, map[1<<15], mapb[1<<15], dup[NDG6] = {0}, mul[NDG6] = {0};
  unsigned lnk, lnkb;
  struct {
    int sc;    /* star content */
    int degen; /* n!/degen is the symmetry number */
    int npr;   /* number of wiggly edges */
    int id[16][2];
    const char name[16];
  } db[NDG6] = {
    /* cf. ``Reformulation of the Virial Series for Classical Fluids'' Fig. 2
     * F. H. Ree and W. G. Hoover, J. Chem. Phys. 41, 6, 1635-1645.
     * The diagrams have been reordered according to Table II. in
     * ``Virial Coefficients for hard spheres and hard disks'',
     * F. H. Ree and W. G. Hoover, J. Chem. Phys. 40, 4 939-950 */
    /*  0 */ {  24,   1,  0, {{0, 0}}, "15"},
    /*  1 */ { -12,  45,  2, {{0, 1}, {2, 3}}, "13"}, /* two lines */

    /*  2 */ {  16,  15,  3, {{0, 1}, {2, 3}, {4, 5}}, "12a"}, /* three lines */
    /*  3 */ {   8, 180,  3, {{0, 1}, {2, 3}, {3, 4}}, "12b"}, /* a line and a 2-line */

    /*  4 */ {  -4,  90,  4, {{0, 1}, {1, 2}, {3, 4}, {4, 5}}, "11a"}, /* two 2-lines */
    /*  5 */ {  -5, 180,  4, {{0, 1}, {2, 3}, {3, 4}, {4, 5}}, "11b"}, /* a line and a 3-line */
    /*  6 */ {  -4,  60,  4, {{0, 1}, {2, 3}, {3, 4}, {4, 2}}, "11c"}, /* a line and a triangle */
    /*  7 */ {  -6,  60,  4, {{0, 1}, {2, 3}, {2, 4}, {2, 5}}, "11d"}, /* a line and an arrow */

    /*  8 */ {  -1, 360,  5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}}, "10a"}, /* a 5-line */
    /*  9 */ {   4,  45,  5, {{0, 1}, {2, 3}, {3, 4}, {4, 5}, {5, 2}}, "10b"}, /* a line and a 4-ring */
    /* 10 */ {  -4,  72,  5, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}, "10c"}, /* a 5-ring */
    /* 11 */ {   3, 180,  5, {{0, 1}, {2, 3}, {3, 4}, {4, 2}, {2, 5}}, "10d"}, /* a line and an umbrella */

    /* 12 */ {   4,  10,  6, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}}, "9a"}, /* two triangles */
    /* 13 */ {   4,  60,  6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}, "9b"}, /* 6-ring */
    /* 14 */ {   1, 360,  6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 3}}, "9c"}, /* a triangle with a 3-tail */
    /* 15 */ {   3, 360,  6, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {0, 5}}, "9d"}, /* 5-ring with a tail */
    /* 16 */ {  -2,  90,  6, {{0, 1}, {2, 3}, {3, 4}, {4, 5}, {5, 2}, {2, 4}}, "9e"}, /* a line and a diamond */

    /* 17 */ {  -2, 360,  7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 2}}, "8a"}, /* a 6-ring with 0-2 */
    /* 18 */ {  -1,  90,  7, {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}}, "8b"}, /* two connected triangles */
    /* 19 */ {  -2, 180,  7, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}, {5, 0}, {5, 2}}, "8c"},

    /* 20 */ {   1,  15,  7, {{0, 1}, {2, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}}, "8d"},
    /* the above diagram does not use the chain topology,
     * 1-3-4-5 is fully-connected by wiggly lines */

    /* 21 */ {   1, 180,  8, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 2}, {1, 3}},"7"},
    /* 22 */ {   1,  60,  9, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 2}, {3, 5}, {1, 4}}, "6"},
  };
  int idchain[][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};
  int idclaw[][2]  = {{0, 1}, {1, 2}, {2, 3}, {2, 4}, {2, 5}};

  for (i = 0; i < (1 << 15); i++)
    map[i] = mapb[i] = -1; /* invalid entry */

  /* compute multiplicity and build topology */
  for (k = 0; k < NDG6; k++) {
    mul[k] = db[k].sc * db[k].degen * ((15 - db[k].npr) % 2 ? -1 : 1);
    dup[k]  = mktop(6, db[k].npr, db[k].id, k, map,  5, idchain, 0);
    dup[k] += mktop(6, db[k].npr, db[k].id, k, mapb, 5, idclaw, 0);
  }

  rvn_zero(v[0]);
  rvn_zero(vb[0]);
  for (t = 1; t <= nsteps; t++) {
    /* spanning tree is a chain */
    for (i = 0; i < 5; i++)
      rvn_add(v[i + 1], v[i], rvn_rndball0(u));
    /* spanning tree is a claw */
    rvn_copy(vb[1], v[1]);
    rvn_copy(vb[2], v[2]);
    rvn_copy(vb[3], v[3]);
    rvn_add(vb[4], vb[2], rvn_rndball0(u));
    rvn_add(vb[5], vb[2], rvn_rndball0(u));

    /* encode the connectivity */
    for (lnk = lnkb = 0, ipr = 0, i = 0; i < 5; i++)
      for (j = i + 1; j < 6; j++) {
        if (rvn_dist2(v[i], v[j])   >= 1) lnk  |= 1 << ipr;
        if (rvn_dist2(vb[i], vb[j]) >= 1) lnkb |= 1 << ipr;
        ipr++;
      }

    if ((k = map[lnk]) >= 0) acc[k]++;  /* chain backbone */
    if ((k = mapb[lnkb]) >= 0) acc[k]++;  /* claw backbone */
  }

  norm1 = -5./720 * pow(2, 5);
  norm2 = -5./720 * pow(8, 5);
  for (tot = 0, k = 0; k < NDG6; k++) {
    x = 1. * acc[k] / nsteps / dup[k];
    tot += y = x * mul[k];
    printf("%2d A6_%-3s: acc %9.7f, dup %2d sc %+3d, sc*dg %+5d mul %+5d "
        "x*mul %+10.7f (%+9.6f,%+10.5f), tot %+10.7f (%+9.6f,%+10.5f)\n",
        k, db[k].name, x, dup[k], db[k].sc, db[k].sc * db[k].degen,
        mul[k], y, norm1 * y, norm2 * y, tot, norm1 * tot, norm2 * tot);
  }
  tot *= -5./720;
  printf("vir6: %.8f, %.7f, %.6f\n", tot, tot * pow(2, 5), tot * pow(8, 5));
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
  printf("D %d, n %d, %ld Monte Carlo moves\n", D, n, nsteps);

  if (n == 3) vir3(nsteps);
  else if (n == 4) vir4(nsteps);
  else if (n == 5) vir5(nsteps);
  else if (n == 6) vir6(nsteps);
  mtsave(NULL);
  return 0;
}

