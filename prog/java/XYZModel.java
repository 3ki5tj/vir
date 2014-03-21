import java.awt.*;



/** A set of atoms with 3D coordinates */
class XYZModel {
  double realXYZ[][]; // 3D real coordinates [np][3]
  int screenXYZ[][];  // 3D screen coordinates in pixels [np][3]
                    // only the first two dimensions are used in drawing
  int zOrder[]; // z-order, larger is closer to the viewer
  int np = -1;

  boolean transformed;
  // rotation/scaling/translation matrix for the conversion from realXYZ[] to screenXYZ[]
  Matrix3D mat = new Matrix3D();
  double real2Screen = 50.0; // real size --> screen size
  double ballSize = 0.5; // ball size (radius) in terms of the real coordinates
                         // 0.5 for hard spheres

  Atom atomDefault = new Atom(0.1, 0.7, 0.1, 1.0, 0.5, 0.5, 1.0);
  Atom atoms[]; // for colors of atoms

  XYZModel() {}

  /** Set the color of a particular atom */
  void setAtom(int i, Atom atom) {
    if ( i >= 0 && i < atoms.length )
      atoms[i] = atom;
  }

  /** Refresh coordinates
   *  x[0..n-1][3] is a three-dimensional vector
   *  n can be less than x.length */
  void updateXYZ(double x[][], int n, boolean center) {
    if (n != np) {
      //System.out.printf("XYZModel.updateXYZ n %d --> %d, (%g, %g, %g) (%g, %g, %g)\n", np, n, x[0][0], x[0][1], x[0][2], x[n-1][0], x[n-1][1], x[n-1][2]);
      np = n;
      realXYZ = new double [np][3];
      screenXYZ = new int [np][3];
      zOrder = new int [np];
      atoms = new Atom [np];
      for (int i = 0; i < np; i++)
        atoms[i] = atomDefault;
    }

    for (int d = 0; d < 3; d++) {
      double xc = 0;
      if ( center ) {
        for (int i = 0; i < np; i++)
          xc += x[i][d];
        xc /= np;
      }
      for (int i = 0; i < np; i++)
        realXYZ[i][d] = x[i][d] - xc;
    }

    transformed = false;
  }

  /** Set the view matrix
   *  `s' is the scaling factor of translating real coordinates
   *    to the screen coordinates
   *  (x0, y0) the screen coordinates of the center */
  void setMatrix(Matrix3D viewMat, double s, double x0, double y0) {
    mat.unit();
    mat.mult(viewMat);
    mat.scale(s, s, s);
    real2Screen = s;
    mat.translate(x0, y0, 0);
    transformed = false;
  }

  /** Get the span of the model
   *  `n' may be less than x.length */
  double getSpan(double x[][], int n) {
    int dim = x[0].length;
    double realSpan = 0, del, fw, fh;

    for (int d = 0; d < dim; d++) {
      double xmin = 1e30, xmax = -1e30;
      for (int i = 0; i < n; i++)
        if ( x[i][d] < xmin ) xmin = x[i][d];
        else if ( x[i][d] > xmax ) xmax = x[i][d];
      if ( (del = xmax - xmin) > realSpan )
        realSpan = del;
    }
    return realSpan;
  }

  /** Translate a given span, return the real-to-screen ratio
   *  `w' and `h' are the width and height of the screen in pixels */
  double getScaleFromSpan(double realSpan, int w, int h) {
    realSpan += ballSize * 2; // add two radii
    double fw = w / realSpan;
    double fh = h / realSpan;
    double facShrink = 0.9; // shrink a bit for the margin
    return (fw < fh ? fw : fh) * facShrink;
  }

  /** Compute the Z-order */
  void getZOrder() {
    // transform the coordinates
    if ( !transformed ) {
      mat.transform(realXYZ, screenXYZ, np);
      transformed = true;
    }

    // bubble sort z-order
    // zOrder[0] is the fartherest from the viewer
    // zOrder[np - 1] is the nearest to the viewer
    for (int i = 0; i < np; i++)
      zOrder[i] = i;
    for (int i = 0; i < np; i++) {
      // find the particle with the smallest z
      int jm = i, k;
      for (int j = i + 1; j < np; j++)
        if (screenXYZ[zOrder[j]][2] < screenXYZ[zOrder[jm]][2])
          jm = j;
      if (jm != i) {
        k = zOrder[i];
        zOrder[i] = zOrder[jm];
        zOrder[jm] = k;
      }
    }
  }

  /** Draw atom */
  void drawAtom(Graphics g, int iz) {
    int i = zOrder[iz];
    Atom atom = atoms[i];
    if (atom == null) return;
    double zMin = screenXYZ[zOrder[0]][2]; // fartherest from the viewer
    double zMax = screenXYZ[zOrder[np - 1]][2]; // closest to the viewer
    int greyScale = Atom.nBalls - 1;
    if (zMin != zMax) // zMin == zMax means the two-dimensional case
      greyScale = (int) (Atom.nBalls * (screenXYZ[i][2] - zMin) / (zMax - zMin) - 1e-6);
    // the atom closest to the viewer has a greyScale of Atom.nBalls - 1
    // the atom fartherest from the viewer has a greyScale of 0
    double radius = ballSize * atom.relRadius * real2Screen;
    atom.paint(g, screenXYZ[i][0], screenXYZ[i][1], greyScale, radius);
  }

  /** Paint this model to the graphics context `g' */
  void paint(Graphics g) {
    if (realXYZ == null || np <= 0) return;
    getZOrder();
    for (int iz = 0; iz < np; iz++)
      drawAtom(g, iz);
  }
}

