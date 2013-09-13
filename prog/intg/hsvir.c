/* first few virial coefficients of the 3D hard spheres system
 * by direct Monte Carlo integration of Mayer diagrams
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
  long t, acc44 = 0, acc45 = 0, acc46 = 0;
  int lnk13, lnk02;
  rvn_t u, v[4];
  double x, rat;

  for (t = 1; t <= nsteps; t++) {
    rvn_zero(v[0]);
    rvn_add(v[1], v[0], rvn_rndball0(u));
    rvn_add(v[2], v[1], rvn_rndball0(u));
    rvn_add(v[3], v[2], rvn_rndball0(u));
    if (rvn_sqr(v[3]) < 1) {
      acc44++;
      lnk02 = (rvn_sqr(v[2]) < 1);
      lnk13 = (rvn_dist2(v[1], v[3]) < 1);
      acc45 += lnk02 + lnk13;
      if (lnk02 && lnk13) acc46++;
    }
  }
  acc45 /= 2; /* for the double counting */
  printf("acc44 %g, acc45 %g, acc46 %g\n", 1.*acc44/nsteps,
      1.*acc45/nsteps, 1.*acc46/nsteps);
  x = -3 * (acc44/8. - acc45/4. + acc46/24.) / nsteps;
  rat = (acc44*3. - acc45*6. + acc46)/(acc44*3. + acc45*6. + acc46);
  printf("vir4: %.8f, %.7f, %.9f, rate %g\n", x, x * pow(2, 3), x * pow(8, 3), rat);
  return x;
}



/* fifth-order virial coefficients */
static double vir5(long nsteps)
{
  long t, acc[10];
  rvn_t u, v[5];
  double tot, abstot, rat, x, y;
  int i, S[10] = {-10, 2, 12, -4, -2,  -12, 4, 8, -12, 120}, lnk[5][5];
  const char *names[10] = {"A55", "A56\'", "A56\'\'", "A57\'", "A57\'\'",
    "A57\'\'\'", "A58\'", "A58\'\'", "A59", "A510"};

  for (i = 0; i < 10; i++) acc[i] = 0;
  rvn_zero(v[0]);
  for (t = 1; t <= nsteps; t++) {
    rvn_add(v[1], v[0], rvn_rndball0(u));
    rvn_add(v[2], v[1], rvn_rndball0(u));
    rvn_add(v[3], v[2], rvn_rndball0(u));
    rvn_add(v[4], v[3], rvn_rndball0(u));

#define COMPLINK(i, j) lnk[i][j] = (rvn_dist2(v[i], v[j]) < 1)
    COMPLINK(0, 2);
    COMPLINK(0, 3);
    COMPLINK(0, 4);
    COMPLINK(1, 3);
    COMPLINK(1, 4);
    COMPLINK(2, 4);

    if (lnk[0][4]) { /* ring diagrams */
      acc[0]++;
      if (lnk[0][3]) {
        acc[1]++;
        if (lnk[2][4]) {
          acc[3]++;
          if (lnk[0][2]) {
            acc[6]++;
            if (lnk[1][3]) {
              acc[8]++;
              if (lnk[1][4]) {
                acc[9]++;
              }
            }
          }
        }
      }
      if (lnk[0][2] && lnk[2][4]) {
        acc[4]++;
        if (lnk[1][3])
          acc[7]++;
      }
    }

    /* the two diamond diagrams */
    if (lnk[1][4] && lnk[0][3]) {
      acc[2]++;
      if (lnk[1][3])
        acc[5]++;
    }
  }
  for (abstot = tot = 0, i = 0; i < 10; i++) {
    double norm1 = -4 * pow(2, 4), norm2 = -4 * pow(8, 4);

    x = 1. * acc[i] / nsteps;
    tot += y = x / S[i];
    abstot += fabs(y);
    printf("%-6s acc%d: x %10.8f, S %3d, x/S %12.8f (%12.7f, %12.6f), "
        "tot %12.8f (%12.7f, %12.6f)\n",
        names[i], i, x, S[i], y, norm1 * y, norm2 * y,
        tot, norm1 * tot, norm2 * tot);
  }
  rat = tot/abstot;
  tot *= -4;
  printf("vir5: %.8f, %.7f, %.6f rate %g\n", tot, tot * pow(2, 4), tot * pow(8, 4), rat);
  return tot;
}



int main(int argc, char *argv[])
{
  long nsteps = 10000000;
  int n = 5;

  /* Metropolis. Rosenbluth, Roksenbluth, Teller & Teller, J. Chem. Phys. 21, 6, 1087-1092, 1953
   * Ree, Hoover, J. Chem. Phys. 46, 4, 939-950, 1964
   * http://www.sklogwiki.org/SklogWiki/index.php/Hard_sphere:_virial_coefficients
   * B2 = b0 = 2 pi/3
   * B3 = 5/8 b0^2
   * B4 = 0.2869495 b0^3
   * B5 = 0.110252(1) b0^4
   * B6 = 0.03888198(91) b0^5
   * B7 = 0.01302354(91) b0^6
   * B8 = 0.0041832(11) b0^7
   * B9 = 0.0013094(13) b0^8
   * B10 = 0.0004035(15) b0^9 */
  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nsteps = strtol(argv[2], NULL, 10);
  printf("n %d, %ld Monte Carlo moves\n", n, nsteps);
  if (n == 3) vir3(nsteps);
  else if (n == 4) vir4(nsteps);
  else if (n == 5) vir5(nsteps);
  return 0;
}

