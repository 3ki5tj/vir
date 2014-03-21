import java.awt.*;



/** XYZModel for a typical MC simulation */
class XYZModelMC extends XYZModel {
  int moveAtom = -1;
  boolean moveAcc = false;

  // failed move
  Atom atomTrial  = new Atom(1.0, 0.1, 0.1, 1.0, 0.5, 0.3, 2.0); // failed trial postion
  Atom atomFailed = new Atom(0.0, 0.0, 0.4, 1.0, 0.5, 0.3, 2.0); // new and old position
  // successful move
  Atom atomMoved  = new Atom(0.1, 1.0, 0.1, 1.0, 0.7, 0.3, 2.0); // trial and new position
  Atom atomOldPos = new Atom(0.2, 0.2, 0.2, 1.0, 0.5, 0.5, 2.0); // old position

  /** Constructor from the real coordinates
   *  x[][3] is a three-dimensional vector */
  XYZModelMC() {
    atomDefault = new Atom(0.1, 0.2, 1.0, 1.0, 0.6, 0.4, 1.0);
  }

  /** Set the atom */
  public void setMoveAtom(int mvAtom, boolean mvAcc, Atom atms[]) {
    moveAtom = mvAtom;
    moveAcc = mvAcc;
    if (mvAtom < 0) return;
    for (int i = 0; i < np; i++)
      atoms[i] = (atms != null) ? atms[i] : atomDefault;
    if ( moveAtom < 0 || moveAtom >= np - 1 ) {
      atoms[np - 1] = null;
    } else {
      atoms[moveAtom] = moveAcc ? atomMoved : atomFailed;
      atoms[np - 1] = moveAcc ? atomOldPos : atomTrial;
    }
  }
}

