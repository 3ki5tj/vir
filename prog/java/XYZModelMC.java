import java.awt.*;



enum AtomTypeMC {
  ATOM_NORMAL, ATOM_TRIAL, ATOM_FAILED, ATOM_MOVED, ATOM_OLDPOS
}

/** XYZModel for a typical MC simulation */
class XYZModelMC extends XYZModel {
  int moveAtom = -1;
  boolean moveAcc = false;

  // normal atom
  Atom atomNormal = new Atom(0.1, 0.2, 1.0, 1.0, 0.5, 0.5, 2.0);
  // failed move
  Atom atomTrial  = new Atom(1.0, 0.2, 0.1, 1.0, 0.5, 0.5, 2.0);  // failed trial postion
  Atom atomFailed = new Atom(0.0, 0.1, 0.3, 1.0, 0.5, 0.5, 2.0);  // new and old position
  // successful move
  Atom atomMoved  = new Atom(0.1, 0.6, 0.7, 1.0, 0.5, 0.5, 2.0);  // trial and new position
  Atom atomOldPos = new Atom(0.2, 0.2, 0.2, 1.0, 0.5, 0.5, 2.0);  // old position

  /** Constructor from the real coordinates
   *  x[][3] is a three-dimensional vector */
  XYZModelMC(double x[][], int n, boolean center) {
    //System.out.printf("XYZModelMC Constructor: n %d, (%g, %g, %g) (%g, %g, %g)\n", n, x[0][0], x[0][1], x[0][2], x[n-1][0], x[n-1][1], x[n-1][2]);
    setAtomDefault(atomNormal);
    updateXYZ(x, n, center);
  }

  /** Set the color of a certain type */
  void setColor(AtomTypeMC colorType, int r, int g, int b) {
    switch ( colorType ) {
      case ATOM_NORMAL:
        atomNormal = new Atom(r, g, b);
        setAtomDefault(atomNormal);
        break;
      case ATOM_TRIAL:  atomTrial  = new Atom(r, g, b); break;
      case ATOM_FAILED: atomFailed = new Atom(r, g, b); break;
      case ATOM_MOVED:  atomMoved  = new Atom(r, g, b); break;
      case ATOM_OLDPOS: atomOldPos = new Atom(r, g, b); break;
      default: break;
    }
  }

  /** Copy coordinates from MCSamp to XYZModel
   *  x[][3] is the wrapped coordinates
   *  x[x.length - 1] is the trial position */
  public void update(double x[][], boolean center,
      int mvAtom, boolean mvAcc) {
    int n = x.length;
    //System.out.printf("XYZModelMC.update np %d, n %d, (%g, %g, %g) (%g, %g, %g)\n", np, n, x[0][0], x[0][1], x[0][2], x[n-1][0], x[n-1][1], x[n-1][2]);
    updateXYZ(x, x.length, center);
    moveAtom = mvAtom;
    moveAcc = mvAcc;
    for (int i = 0; i < np; i++)
      atoms[i] = atomNormal;
    if ( moveAtom < 0 || moveAtom >= np - 1 ) {
      atoms[np - 1] = null;
    } else {
      atoms[moveAtom] = moveAcc ? atomMoved : atomFailed;
      if ( moveAcc ) { // the original position before pos
        atoms[np - 1] = atomOldPos;
      } else { // the failed trial position
        atoms[np - 1] = atomTrial;
      }
    }
  }
}

