import java.awt.*;



/** A set of atoms with 3D coordinates */
class XYZModel {
  double realXYZ[][]; // 3D real coordinates [np][3]
  int screenXYZ[][];  // 3D screen coordinates in pixels [np][3]
                    // only the first two dimensions are used in drawing
  int zorder[]; // z-order, larger is closer to the viewer
  int np = -1;

  boolean transformed;
  Matrix3D mat; // rotation/scaling/translation matrix for the conversion
                // from realXYZ[] to screenXYZ[]
  double real2Screen = 100.0; // real size --> screen size
  double ballSize = 0.5; // ball size (radius) in terms of the real coordinates
                          // it is 0.5 because we are simulating hard spheres
  int moveAtom = -1;
  boolean moveAcc = false;

  // various atom types
  static Atom atomNormal = new Atom(  0,   0, 255);
  // failed move
  static Atom atomTrial  = new Atom(255,   0,   0);  // failed trial postion
  static Atom atomFailed = new Atom( 50,   0, 100);  // new and old position
  // successful move
  static Atom atomMoved  = new Atom(  0, 255,   0);  // trial and new position
  static Atom atomOldpos = new Atom(  0, 100, 100);  // old position

  /** Constructor from the real coordinates
   *  r[][3] is a three-dimensional vector */
  XYZModel(double r[][], boolean centering) {
    mat = new Matrix3D();
    updateXYZ(r, centering);
  }

  /** Refresh coordinates
   *  r[][3] is a three-dimensional vector */
  void updateXYZ(double r[][], boolean centering) {
    int n = r.length;
    if (n != np) {
      np = n;
      realXYZ = new double[np][3];
      screenXYZ = new int[np][3];
      zorder = new int[np];
    }

    for (int d = 0; d < 3; d++) {
      double xc = 0;
      if ( centering ) {
        for (int i = 0; i < np; i++)
          xc += r[i][d];
        xc /= np;
      }
      for (int i = 0; i < np; i++)
        realXYZ[i][d] = r[i][d] - xc;
    }

    transformed = false;
  }

  /** Set of the view matrix
   *  s is the scaling factor of translating real coordinates
   *    to the screen coordinates
   *  (x0, y0) the screen coordinates of the center */
  void setMatrix(Matrix3D amat, double s, double x0, double y0) {
    mat.unit();
    mat.mult(amat);
    mat.scale(s, s, s);
    real2Screen = s;
    mat.translate(x0, y0, 0);
    transformed = false;
  }

  /** Paint this model to the graphics context `g' */
  void paint(Graphics g) {
    if (realXYZ == null || np <= 0) return;

    // transform the coordinates
    if ( !transformed ) {
      mat.transform(realXYZ, screenXYZ, np);
      transformed = true;
    }

    // bubble sort z-order
    // zorder[0] is the fartherest from the viewer
    // zorder[np - 1] is the nearest to the viewer
    for (int i = 0; i < np; i++)
      zorder[i] = i;
    for (int i = 0; i < np; i++) {
      // find the particle with the smallest z
      int jm = i, k;
      for (int j = i + 1; j < np; j++)
        if (screenXYZ[zorder[j]][2] < screenXYZ[zorder[jm]][2])
          jm = j;
      if (jm != i) {
        k = zorder[i];
        zorder[i] = zorder[jm];
        zorder[jm] = k;
      }
    }

    double zmin = screenXYZ[zorder[0]][2];
    double zmax = screenXYZ[zorder[np - 1]][2];
    //System.out.println("" + zmin + " " + zmax);

    // radius in pixels
    double radius = real2Screen * ballSize;

    for (int iz = 0; iz < np; iz++) {
      int i = zorder[iz], greyScale;
      if (zmin == zmax) { // two-dimensional case
        greyScale = 15;
      } else {
        greyScale = (int) (15.999 * (screenXYZ[i][2] - zmin) / (zmax - zmin));
      }
      // screenXYZ[3*i] and screenXYZ[3*i + 1] are the (x, y) coordinates on screen
      Atom atom = atomNormal;
      if ( i == moveAtom )
        atom = moveAcc ? atomMoved : atomFailed;
      if ( i == np - 1 ) {
        if ( moveAtom < 0 ) {
          atom = null;
        } else if ( moveAcc ) { // the original position before pos
          atom = atomOldpos;
        } else { // the failed trial position
          atom = atomTrial;
        }
      }
      if ( atom != null )
        atom.paint(g, screenXYZ[i][0], screenXYZ[i][1],
                   greyScale, radius);
    }
  }
}




