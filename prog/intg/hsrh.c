/* first few virial coefficients of the 3D hard spheres system
 * by direct Monte Carlo integration of Ree-Hoover diagrams
 * Copyright (C) 2013 Cheng Zhang */
#include <time.h>



/* #define `D' here to simulate dimensions other than D = 3 */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"


#ifdef MPI
#include <mpi.h>
MPI_Comm comm = MPI_COMM_WORLD;
#endif
#define MASTER 0
int inode = MASTER, nnodes = 1;


long nstrep = 100000000;


/* make an output file name */
#define mkfnout(fn, n) { sprintf(fn, "intD%dn%d.dat%d", D, n, inode); \
  if (inode == MASTER) fn[strlen(fn) - 1] = '\0'; }



/* compute the virial coefficient and save it to file */
static double save(int n, long nsteps, int m, const long *acc,
    const double *fac, const char **names)
{
  FILE *fp;
  int i;
  double vir = 0;
  char fnout[64] = "";

  for (i = 0; i < m; i++)
    vir += acc[i] * fac[i];
  vir /= nsteps;
  for (i = 1; i <= n; i++) vir /= i;
  vir *= -(n - 1) * pow(2, n - 1);

  mkfnout(fnout, n);
  xfopen(fp, fnout, "w", return vir);
  fprintf(fp, "#0 %d %d %d V0\n%ld\n", D, n, m, nsteps);
  for (i = 0; i < m; i++)
    fprintf(fp, "%16ld %+.14e %s\n", acc[i], fac[i], names[i]);
  fprintf(fp, "%.14e\n", vir);
  fclose(fp);
  fprintf(stderr, "%4d: saving file to %s\n", inode, fnout);
  if (inode == MASTER) mtsave(NULL);
  return vir;
}



/* third-order virial coefficients */
static double vir3(long nsteps)
{
  long t, acc = 0;
  rvn_t v, u;
  double x;

  rvn_zero(u);
  for (t = 1; t <= nsteps; t++) {
    u[0] = pow(rnd0(), 1./D);
    rvn_rndball0(v); v[0] += u[0];
    if (rvn_sqr(v) < 1) acc++;
  }
  x = -2 * (-acc/6.) / nsteps;
  printf("vir3: %.9e, %.9e, %.9e, acc %g\n",
      x, x * pow(2, 2), x * pow(8, 2), 1. * acc / nsteps);
  return x;
}



/* fourth-order virial coefficients */
static double vir4(long nsteps)
{
  long t, acc[2] = {0, 0};
  int lnk02, lnk13;
  double fac[2] = {3, -2}, vir = 0;
  const char *names[] = {"A4_4", "A4_6"};
  rvn_t u, v[4];

  rvn_zero(v[0]);
  rvn_zero(v[1]);
  for (t = 1; t <= nsteps; t++) {
    v[1][0] = pow(rnd0(), 1./D);
    rvn_rndball0(v[2]); v[2][0] += v[1][0];
    rvn_add(v[3], v[2], rvn_rndball0(u));
    if (rvn_sqr(v[3]) < 1) {
      lnk02 = ( rvn_sqr(v[2]) < 1 );
      lnk13 = ( rvn_dist2(v[1], v[3]) < 1 );
      if (lnk02 && lnk13) acc[1]++;
      else if (!lnk02 && !lnk13) acc[0]++;
    }
    if (t % nstrep == 0 || t == nsteps)
      vir = save(4, t, 2, acc, fac, names);
  }
  printf("B4/B2^3: %.9e\n", vir);
  return vir;
}



#define NDG5 5

/* fifth-order virial coefficients */
static double vir5(long nsteps)
{
  rvn_t u, v[5];
  double fac[NDG5], vir = 0;
  long t, acc[NDG5 + 1] = {0}; /* acc[NDG5] is reserved for invalid diagrams */
  int sc[NDG5] = {-6, 3, -2, 1, 1}, degen[NDG5] = {1, 15, -30, 10, -12};
  int i, mul[NDG5], dup[NDG5] = {1, 6, 7, 1, 1}, map[64];
  const char *names[NDG5] = {"A5_10", "A5_8", "A5_7", "A5_6", "A5_5"};
  unsigned code;

  for (i = 0; i < NDG5; i++) {
    mul[i] = sc[i] * degen[i];
    fac[i] = 1. * mul[i] / dup[i];
  }

  for (i = 0; i < 64; i++)
    map[i] = NDG5; /* `NDG5' means invalid or un-biconnected diagrams */
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
  rvn_zero(v[1]);
  for (t = 1; t <= nsteps; t++) {
    v[1][0] = pow(rnd0(), 1./D);
    rvn_rndball0(v[2]); v[2][0] += v[1][0];
    rvn_add(v[3], v[2], rvn_rndball0(u));
    rvn_add(v[4], v[3], rvn_rndball0(u));

#define CNT(i, j) (rvn_dist2(v[i], v[j]) < 1)
    code = CNT(0, 4)       | (CNT(0, 2) << 1) | (CNT(0, 3) << 2)
        | (CNT(1, 3) << 3) | (CNT(1, 4) << 4) | (CNT(2, 4) << 5);

    acc[ map[code] ]++;
    if (t % nstrep == 0 || t == nsteps)
      vir = save(5, t, NDG5, acc, fac, names);
  }
  printf("B5/B2^4: %.9e\n", vir);
  return vir;
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
  double fac[NDG6], vir = 0;
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
  const char *names[NDG6];

  for (i = 0; i < (1 << 15); i++)
    map[i] = mapb[i] = -1; /* invalid entry */

  /* compute multiplicity and build topology */
  for (k = 0; k < NDG6; k++) {
    names[k] = db[k].name;
    mul[k] = db[k].sc * db[k].degen * ((15 - db[k].npr) % 2 ? -1 : 1);
    dup[k]  = mktop(6, db[k].npr, db[k].id, k, map,  5, idchain, 0);
    dup[k] += mktop(6, db[k].npr, db[k].id, k, mapb, 5, idclaw, 0);
    fac[k] = 1.*mul[k]/dup[k];
  }

  rvn_zero(v[0]);
  rvn_zero(v[1]);
  rvn_zero(vb[0]);
  rvn_zero(vb[1]);
  for (t = 1; t <= nsteps; t++) {
    /* spanning tree is a chain */
    v[1][0] = pow(rnd0(), 1./D);
    rvn_rndball0(v[2]); v[2][0] += v[1][0];
    rvn_add(v[3], v[2], rvn_rndball0(u));
    rvn_add(v[4], v[3], rvn_rndball0(u));
    rvn_add(v[5], v[4], rvn_rndball0(u));
    /* spanning tree is a claw */
    vb[1][0] = v[1][0];
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
    if (t % nstrep == 0 || t == nsteps)
      vir = save(6, t, NDG6, acc, fac, names);
  }
  printf("B6/B2^5: %.14e\n", vir);
  return vir;
}



int main(int argc, char *argv[])
{
  long nsteps = 10000000;
  int n = 6;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &inode);
  MPI_Comm_size(comm, &nnodes);
#endif
  /* Metropolis, Rosenbluth, Roksenbluth, Teller & Teller, J. Chem. Phys. 21, 6, 1087-1092, 1953
   * Ree, Hoover, J. Chem. Phys. 46, 4, 939-950, 1964
   * http://www.sklogwiki.org/SklogWiki/index.php/Hard_sphere:_virial_coefficients
   * B2  = b0 = 2 pi/3
   * B3  = 5/8 b0^2
   * B4  = 0.2869495... b0^3
   * B5  = 0.11025175(4) b0^4
   * B6  = 0.0388823(3)  b0^5
   * B7  = 0.0130228(4)  b0^6
   * B8  = 0.0041832(11) b0^7
   * B9  = 0.0013092(12) b0^8
   * B10 = 0.0004035(15) b0^9 */
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nsteps = (long) (atof(argv[2]) + .5);
  if (argc > 3) nstrep = (long) (atof(argv[3]) + .5);
  mtscramble(inode * 2034091783u + time(NULL));
  if (inode == MASTER)
    printf("D %d, n %d, %ld Monte Carlo moves\n", D, n, nsteps);

  if (n == 3) vir3(nsteps);
  else if (n == 4) vir4(nsteps);
  else if (n == 5) vir5(nsteps);
  else if (n == 6) vir6(nsteps);
#ifdef MPI
  MPI_Finalize();
#endif
  return 0;
}

