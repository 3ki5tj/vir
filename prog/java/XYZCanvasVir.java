import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;



/** Panel for drawing particles */
class XYZCanvasVir extends XYZCanvas {
  public static final long serialVersionUID = 5L;

  boolean useMDS = false; // use multidimensional scaling to transform coordinates
  boolean colorMDS = false; // use MDS to color atoms
  Atom atoms[];

  /** Get 3D position
   *  `n' may be less than x.length */
  private double [][] getPos3D(double x[][], int n, double xMove[]) {
    int D = x[0].length;
    double xyz[][] = new double[n + 1][3];

    // copy the direct coordinates
    for (int i = 0; i < n; i++) {
      xyz[i][0] = x[i][0];
      xyz[i][1] = x[i][1];
      xyz[i][2] = (D == 2) ? 0.0 : x[i][2];
    }
    if (xMove == null) xMove = x[0];
    // the last is the trial position
    xyz[n][0] = xMove[0];
    xyz[n][1] = xMove[1];
    xyz[n][2] = (D == 2) ? 0 : xMove[2];

    if ( useMDS ) {
      double x0[][] = new double[n + 1][D];

      // copy coordinates with the additional trial position
      for (int i = 0; i < n; i++)
        for (int k = 0; k < D; k++)
          x0[i][k] = x[i][k];
      for (int k = 0; k < D; k++)
        x0[n][k] = xMove[k];

      MDS mds = new MDS(x0);
      double eneMin = mds.min(xyz);
      // set the colors of the atoms by the MDS discrepency
      if ( atoms == null || atoms.length < n + 1 )
        atoms = new Atom [n + 1];
      if ( colorMDS ) mds.setAtoms(xyz, atoms);
      //System.out.printf("D0 %d, mds minimal energy: %g\n", D, eneMin);
    }
    return xyz;
  }

  /** Refresh the coordinates
   *  x[][] is the wrapped coordinates
   *  `n' may be less than x.length */
  public void refresh(double x[][], int n, double xo[], double xi[],
      boolean center, int mvAtom, boolean mvAcc, boolean adjScale) {
    if (model == null) {
      model = new XYZModelMC();
      adjScale = true;
    }
    XYZModelMC modelMC = (XYZModelMC) model;
    double xyz[][] = getPos3D(x, n,
        mvAtom < 0 ? null : mvAcc ? xo : xi);
    modelMC.updateXYZ(xyz, n + 1, center);
    modelMC.setMoveAtom(mvAtom, mvAcc,
        (useMDS && colorMDS) ? atoms : null);
    //System.out.printf("%s\n", getSize());
    if ( adjScale ) {
      // empirical formulas for the span of the cluster
      int dim = x[0].length;
      double idealSpan;
      if ( useMDS )
        idealSpan = 1.5 * Math.pow(n, 1./3);
      else
        idealSpan = 1.5 * Math.pow(n*n/40., 1./dim);
      realSpan = max(model.getSpan(x, n), idealSpan);
    }
    repaint();
  }
}

