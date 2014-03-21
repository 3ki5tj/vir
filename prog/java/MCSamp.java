import static java.lang.Math.*;
import java.util.Random;



/** Metropolis sampling */
class MCSamp {
  public int D = 3; // the dimension
  public int N = 7; // the order
  public Diagram g, ng;

  int step; // the current simulation step
  double mcAmp = 1.5; // amplitude of randomly moving a particle
  int sampFreq = 1; // sampling frequency
  long mcAcc = 0, mcTot = 0;
  double x[][], xi[], xo[];
  Ave avFb = new Ave(), avNr = new Ave();
  double xyz[]; // 3D coordinates to output

  double fb, nr; // star and ring contents
  private boolean fbDirty = true;

  Random rng = new Random();

  /** Constructor */
  MCSamp(int dim, int order) {
    init(dim, order);
  }

  /** Initialize system */
  public void init(int dim, int order) {
    D = dim;
    N = order;
    mcAmp = 1.5/D;
    x = new double [N][D];
    xyz = new double [(N + 1) * 3];
    double rr = 0.7/(2.*sin(PI/N));
    for (int i = 0; i < N; i++) {
      x[i][0] = rr * cos(2*i*PI/N);
      x[i][1] = rr * sin(2*i*PI/N);
      for (int j = 2; j < D; j++) x[i][j] = 0;
    }
    xi = new double [D];
    xo = new double [D];
    g  = new Diagram(N, D, x);
    ng = new Diagram(N, D, x);
    clearData();
  }

  /** Clear trajectory data */
  void clearData() {
    step = 0;
    mcAcc = mcTot = 0;
    avFb.clear();
    avNr.clear();
  }

  public int moveAtom = -1;
  public boolean moveAcc = false;

  /** Step of the Metropolis algorithm */
  public void mcStep() {
    int i, j, k;

    moveAtom = i = (int) (N * rng.nextDouble());
    // displace particle i
    for (k = 0; k < D; k++) {
      xo[k] = x[i][k];
      xi[k] = xo[k] + (2 * rng.nextDouble() - 1) * mcAmp;
    }
    ng.copy(g);
    for (j = 0; j < N; j++) {
      if (j == i) continue;
      // distance between i and j
      double dx, x2 = 0;
      for (k = 0; k < D; k++) {
        dx = xi[k] - x[j][k];
        x2 += dx * dx;
      }
      if (x2 < 1) {
        ng.link(i, j);
      } else {
        ng.unlink(i, j);
      }
    }
    mcTot += 1;
    if ( ng.biconnected() ) {
      mcAcc += 1;
      for (k = 0; k < D; k++)
        x[i][k] = xi[k];
      g.copy(ng);
      //for (int ii = 0; ii < N; ii++) System.out.print("x" + ii + ": " + x[ii][0] + " " + x[ii][1] + " " + x[ii][2] + "; ");
      //System.out.println("moving " + i);
      moveAcc = true;
      fbDirty = true;
    } else {
      moveAcc = false;
    }
    if (sampFreq > 0 && (step + 1) % sampFreq == 0) {
      if ( fbDirty ) {
        if (N <= DiagramMap.NMAX) {
          double [] arr2 = g.getFbNrMap();
          fb = arr2[0];
          nr = arr2[1];
        } else {
          fb = g.getFb();
          nr = g.getNr();
        }
        fbDirty = false;
      }
      avFb.add(fb);
      avNr.add(nr);
    }
    step++;
  }

  /** Remove the center of mass motion */
  void rmcom() {
    for (int k = 0; k < D; k++) {
      double xc = 0;
      for (int i = 0; i < N; i++)
        xc += x[i][k];
      xc /= N;
      for (int i = 0; i < N; i++)
        x[i][k] -= xc;
      xi[k] -= xc;
      xo[k] -= xc;
    }
  }

  /** Get the virial coefficients */
  public double getVir() {
    double fb = avFb.getAve();
    double sc = avNr.getAve();
    if (sc == 0) return 0;
    double Bring = abs( Diagram.BringArr[D][N] );
    return -fb/sc*Bring;
  }
}



