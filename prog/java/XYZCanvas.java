import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;



/** Panel for drawing particles */
class XYZCanvas extends JPanel
    implements MouseListener, MouseMotionListener, MouseWheelListener {
  Image img;
  Graphics imgG;
  Dimension imgSize;

  double real2Screen;
  double zoomScale = 1.0;
  public Matrix3D viewMatrix = new Matrix3D(); // view matrix
  private Matrix3D tmpMatrix = new Matrix3D(); // temporary matrix
  int mouseX, mouseY; // mouse position
  String message;
  XYZModelMC model;

  public static final long serialVersionUID = 2L;

  public XYZCanvas() {
    super();
  }

  boolean useMDS = false; // use multidimensional scaling to transform coordinates

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
    // the last is the trial position
    xyz[n][0] = xMove[0];
    xyz[n][1] = xMove[1];
    xyz[n][2] = (D == 2) ? 0 : xMove[2];

    if ( D > 3 && useMDS ) {
      double x0[][] = new double[n + 1][D];

      // copy coordinates with the additional trial position
      for (int i = 0; i < n; i++)
        for (int k = 0; k < D; k++)
          x0[i][k] = x[i][k];
      for (int k = 0; k < D; k++)
        x0[n][k] = xMove[k];

      MDS mds = new MDS(x0);
      double eneMin = mds.min(xyz);
      //System.out.printf("D0 %d, mds minimal energy: %g\n", D, eneMin);
    }
    return xyz;
  }

  /** Copy coordinates from MCSamp to XYZModel
   *  `x' and `xMove' are high dimensional coordinates */
  public void update(double x[][], int n, double xMove[],
                     int moveAtom, boolean moveAcc) {
    model.update(getPos3D(x, n, xMove), false, moveAtom, moveAcc);
  }

  /** Determine scaling factors based on the current coordinates */
  void gauge(double x[][], int n) {
    int dim = useMDS ? 3 : x[0].length;
    model = new XYZModelMC(getPos3D(x, n, x[0]), n + 1, false);
    // empirical formula for the span of the cluster
    double realSpan = 1.5*Math.pow(x.length, 1./dim);
    real2Screen = model.getScaleFromSpan(realSpan, getSize().width, getSize().height);
    if ( newImgBuf() ) repaint();
  }

  public boolean bStarted = false;

  private void start() {
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
    bStarted = true;
  }

  /** Create an image buffer for drawing
   *  return if the canvas is ready */
  public boolean newImgBuf() {
    if (getSize().width == 0 || getSize().height == 0) {
      System.out.println("canvas not ready.");
      return false;
    }
    img = createImage(getSize().width, getSize().height);
    if (imgG != null) { imgG.dispose(); }
    imgG = img.getGraphics();
    imgSize = getSize();
    if ( !bStarted ) start();
    return true;
  }

  public void update(Graphics g) {
    if (img == null)
      g.clearRect(0, 0, getSize().width, getSize().height);
    paintComponent(g);
  }

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    if (model != null) {
      model.setMatrix(viewMatrix, real2Screen * zoomScale,
                      getSize().width/2, getSize().height/2);
      if ( img != null ) {
        if ( !imgSize.equals(getSize()) )
          newImgBuf();
        imgG.setColor(Color.BLACK);
        imgG.fillRect(0, 0, getSize().width, getSize().height);
        model.paint(imgG);
        g.drawImage(img, 0, 0, this);
      } else
        model.paint(g);
    }
  }

  /** Event handling */
  public void mouseClicked(MouseEvent e) { }
  public void mousePressed(MouseEvent e) {
    mouseX = e.getX();
    mouseY = e.getY();
    // consume this event so that it will not be processed
    // in the default manner
    e.consume();
  }
  public void mouseReleased(MouseEvent e) { }
  public void mouseEntered(MouseEvent e) { }
  public void mouseExited(MouseEvent e) { }

  public void mouseDragged(MouseEvent e) {
    int x = e.getX();
    int y = e.getY();
    tmpMatrix.unit();
    tmpMatrix.xrot( 360.0 * (mouseY - y) / getSize().height );
    tmpMatrix.yrot( 360.0 * (x - mouseX) / getSize().width );
    viewMatrix.mult(tmpMatrix);
    repaint();
    mouseX = x;
    mouseY = y;
    e.consume();
  }

  public void mouseMoved(MouseEvent e) { }

  public void mouseWheelMoved(MouseWheelEvent e) {
    int notches = e.getWheelRotation();
    if ((zoomScale -= 0.05f*notches) < 0.09999f)
      zoomScale = 0.1f;
    repaint();
  }
}


