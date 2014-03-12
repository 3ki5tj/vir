import java.awt.*;


/** A set of atoms with 3D coordinates */
class XYZModel {
  double realxyz[][]; // 3D real coordinates 3*np
  int scrnxyz[][];  // 3D screen coordinates in pixels
                    // only the first two dimensions are used in drawing
  int zorder[]; // z-order, larger is closer to the viewer
  int np = -1;

  boolean transformed;
  Matrix3D mat; // rotation/scaling/translation matrix for the conversion
                // from realxyz[] to scrnxyz[]
  double real2scrn = 100.0; // real size --> screen size, to be set by MyCanvas
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

  /** Constructor from the real coordinates */
  XYZModel(double r[][]) {
    mat = new Matrix3D();
    updateXYZ(r);
  }

  /** Refresh coordinates */
  void updateXYZ(double r[][]) {
    int n = r.length;
    if (n != np) {
      np = n;
      realxyz = new double[np][3];
      scrnxyz = new int[np][3];
      zorder = new int[np];
    }
    for (int i = 0; i < np; i++) {
      realxyz[i][0] = r[i][0];
      realxyz[i][1] = r[i][1];
      realxyz[i][2] = r[i][2];
    }
    transformed = false;
    //for (int i = 0; i < np; i++)
    //  System.out.println("i " + i + ": " + r[i][0] + " " + r[i][1] + " " + r[i][2]);
  }

  /** Paint this model to the graphics context `g' */
  void paint(Graphics g) {
    if (realxyz == null || np <= 0) return;

    // transform the coordinates
    if ( !transformed ) {
      mat.transform(realxyz, scrnxyz, np);
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
        if (scrnxyz[zorder[j]][2] < scrnxyz[zorder[jm]][2])
          jm = j;
      if (jm != i) {
        k = zorder[i];
        zorder[i] = zorder[jm];
        zorder[jm] = k;
      }
    }

    double zmin = scrnxyz[zorder[0]][2];
    double zmax = scrnxyz[zorder[np - 1]][2];
    //System.out.println("" + zmin + " " + zmax);

    // radius in pixels
    double radius = real2scrn * ballSize;

    //System.out.println("real2scrn " + real2scrn + ", ball " + ballSize);
    for (int iz = 0; iz < np; iz++) {
      int i = zorder[iz], greyscale;
      if (zmin == zmax) { // two-dimensional case
        greyscale = 15;
      } else {
        greyscale = (int) (15.999 * (scrnxyz[i][2] - zmin) / (zmax - zmin));
      }
      // scrnxyz[3*i] and scrnxyz[3*i + 1] are the (x, y) coordinates on screen
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
        atom.paint(g, scrnxyz[i][0], scrnxyz[i][1], greyscale, radius);
    }
  }
}




