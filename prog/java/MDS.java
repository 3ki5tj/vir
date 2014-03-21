import static java.lang.Math.*;
import java.util.Random;



/** Multidimensional scaling */
class MDS {
  int N = -1;
  double [][] mat;

  /** Constructor from the coordinates */
  MDS(double x[][]) {
    N = x.length;
    mat = new double[N][N];
    getDisMat(mat, x);
  }

  /** Compute the distance matrix */
  void getDisMat(double mat[][], double x[][]) {
    int D = x[0].length; // D can be greater than 3

    for (int i = 0; i < N - 1; i++) {
      for (int j = i + 1; j < N; j++) {
        double d2 = 0;
        for (int k = 0; k < D; k++) {
          double dx = x[i][k] - x[j][k];
          d2 += dx * dx;
        }
        mat[i][j] = mat[j][i] = sqrt(d2);
      }
    }
  }

  int iterMax = 10000;
  double deneTol = 3e-5;
  double eneTol = 0.1;

  /** Get best coordinates */
  double min0(double y[][], double f[][], double yp[][], double fp[][]) {
    if (y.length != N) {
      System.out.printf("dimension mismatch %d vs %d\n", y.length, N);
      System.exit(1);
    }
    int D = y[0].length;
    int npr = N * (N - 1) / 2;

    double ene = getForce(y, f);
    double dt = 0.1;
    int iter;

    for (iter = 0; iter < iterMax; iter++) {
      double enep = ene;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
          // backup the force
          yp[i][j] = y[i][j];
          fp[i][j] = f[i][j];

          // move along the force
          y[i][j] += f[i][j] * dt;
        }
      }
      ene = getForce(y, f);
      //System.out.printf("iter %d, ene %g --> %g\n", iter, enep, ene);
      if (ene > enep) {
        dt *= 0.7;
        // recover the coordinates and force
        for (int i = 0; i < N; i++) {
          for (int j = 0; j < D; j++) {
            y[i][j] = yp[i][j];
            f[i][j] = fp[i][j];
          }
        }
        ene = enep;
      } else {
        if (abs(ene - enep) < deneTol * npr * dt || ene < eneTol)
          break;
        dt *= 1.03f; // attempt to increase the step size
      }
    }
    if (iter >= iterMax)
      System.out.printf("MDS failed in %d steps\n", iter);
    rmcom(y);
    //print(y);
    return ene;
  }

  int numTrials = 10; // number of different starting positions

  /** Get best coordinates */
  double min(double ymin[][]) {
    int D = ymin[0].length;
    double ene0, ene, eneMin = 1e30;

    double [][] y = new double [N][D];
    double [][] f  = new double [N][D];
    double [][] yp = new double [N][D];
    double [][] fp = new double [N][D];

    vecCopy(y, ymin);
    for (int k = 1; k <= numTrials; k++) {
      ene0 = getForce(y, f);
      ene = min0(y, f, yp, fp);
      if (ene < eneMin) {
        eneMin = ene;
        vecCopy(ymin, y);
      }
      //System.out.printf("D %d, trial %2d: ene %g to %g, min %g\n",
      //    y[0].length, k, ene0, ene, eneMin);
      if (eneMin < eneTol) break;
      if (k < numTrials) getRandomPos(y);
    }
    y = yp = f = fp = null;
    return eneMin;
  }

  /** Copy vectors */
  void vecCopy(double a[][], double b[][]) {
    int n = a.length, D = a[0].length, i, j;
    for (i = 0; i < N; i++)
      for (j = 0; j < D; j++)
        a[i][j] = b[i][j];
  }

  Random rng = new Random();

  /** Get a random position */
  void getRandomPos(double y[][]) {
    int D = y[0].length;
    double span = 0;

    for (int j = 0; j < D; j++) {
      double smx = 0, smx2 = 0;
      for (int i = 0; i < N; i++) {
        double yy = y[i][j];
        smx += yy;
        smx2 += yy * yy;
      }
      smx /= N; // average
      smx2 = (smx2 / N - smx * smx); // variance
      span += sqrt(smx2); // standard deviation
    }
    span /= D;

    // assign a random gaussian configuration
    for (int i = 0; i < N; i++)
      for (int j = 0; j < D; j++)
        y[i][j] = span * rng.nextGaussian();
  }

  final double DMIN = 1e-6;

  double diff[];

  /** Compute the force and energy */
  double getForce(double y[][], double f[][]) {
    int D = y[0].length;

    if (diff == null || diff.length != N)
      diff = new double [N];
    for (int i = 0; i < N; i++) diff[i] = 0;

    // clear the force
    if ( f != null )
      for (int i = 0; i < N; i++)
        for (int k = 0; k < D; k++)
          f[i][k] = 0;

    double ene = 0;
    double [] dy = new double [D];
    // compute the pair difference
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double dy2 = 0;
        for (int k = 0; k < D; k++) {
          dy[k] = y[i][k] - y[j][k];
          dy2 += dy[k] * dy[k];
        }
        double dij = sqrt(dy2), dij1 = max(dij, DMIN);
        double dref = mat[i][j], dref1 = max(dref, DMIN);
        double dsc = (dij - dref)/dref1;
        diff[i] += abs(dsc);
        diff[j] += abs(dsc);
        ene += dsc * dsc;
        if (f != null)
          for (int k = 0; k < D; k++) {
            dy[k] *= -(dij - dref) / (dref1 * dref1 * dij1);
            f[i][k] += dy[k];
            f[j][k] -= dy[k];
          }
      }
    }
    dy = null;
    return ene;
  }

  /** Compute the relative discrepency at each atom */
  private double getRelDiff(double x[][]) {
    getForce(x, null);
    double diffmin = 1e30, diffmax = 0;
    for (int i = 0; i < N; i++) {
      if ( diff[i] < diffmin )
        diffmin = diff[i];
      else if ( diff[i] > diffmax )
        diffmax = diff[i];
    }
    if (diffmax > diffmin + 1e-6) {
      for (int i = 0; i < N; i++)
        diff[i] = (diff[i] - diffmin) / (diffmax - diffmin);
    } else {
      for (int i = 0; i < N; i++)
        diff[i] = 0;
    }
    return diffmax;
  }

  /** Set the color according to the relative discrepency */
  void setAtoms(double x[][], Atom atoms[]) {
    getRelDiff(x);
    for (int i = 0; i < N; i++) {
      double w = diff[i], w1, cmin = 0.3;
      if (w < 0.5) {
        w1 = w * 2;
        atoms[i] = new Atom(cmin, cmin + (1 - cmin)*w1, 1 - (1 - cmin)*w1,
            1.0, 0.4, 0.7, 2.0);
      } else {
        w1 = w * 2 - 1;
        atoms[i] = new Atom(cmin + (1 - cmin)*w1, 1 - (1 - cmin)*w1, cmin,
            1.0, 0.4, 0.7, 2.0);
      }
    }
  }

  /** Remove the center of mass */
  void rmcom(double x[][]) {
    int D = x[0].length;
    int N = x.length;

    for (int k = 0; k < D; k++) {
      double xc = 0;
      for (int i = 0; i < N; i++)
        xc += x[i][k];
      xc /= N;
      for (int i = 0; i < N; i++)
        x[i][k] -= xc;
    }
  }

  /** Print coordinates */
  void print(double x[][]) {
    int D = x[0].length;
    int N = x.length;

    for (int i = 0; i < N; i++) {
      System.out.printf("%2d: ", i);
      for (int j = 0; j < D; j++)
        System.out.printf("%+8.3f ", x[i][j]);
      System.out.println("");
    }
  }
}

