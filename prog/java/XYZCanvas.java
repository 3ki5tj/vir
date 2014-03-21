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

  double realSpan; // size of a real span (saved as a copy)
  double zoomScale = 1.0;
  public Matrix3D viewMatrix = new Matrix3D(); // view matrix
  private Matrix3D tmpMatrix = new Matrix3D(); // temporary matrix
  int mouseX, mouseY; // mouse position
  XYZModel model;

  public static final long serialVersionUID = 2L;

  public XYZCanvas() {
    super();
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
  }

  /** Prepare a buffer for the image, return if the canvas is ready
   *  Only need to call this when the size of the canvas is changed,
   *    Since we automatically detect the size change in paint()
   *    only call this on start */
  public boolean newImgBuf() {
    Dimension sz = getSize();
    if (sz.width == 0 || sz.height == 0)
      return false;
    // quit if the current image already has the right size
    if (img != null && imgG != null && sz.equals(imgSize))
      return true;
    img = createImage(sz.width, sz.height);
    if (imgG != null) imgG.dispose();
    imgG = img.getGraphics();
    imgSize = sz;
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
      newImgBuf(); // refresh the image buffer if necessary
      // compute the real-to-screen ratio, this variable differs
      // from model.real2Screen by zoomScale
      Dimension sz = getSize();
      double real2Screen0 = model.getScaleFromSpan(realSpan, sz.width, sz.height);
      model.setMatrix(viewMatrix, real2Screen0 * zoomScale,
                      sz.width/2, sz.height/2);
      imgG.setColor(Color.BLACK);
      imgG.fillRect(0, 0, sz.width, sz.height);
      model.paint(imgG);
      g.drawImage(img, 0, 0, this);
    }
  }

  /** Refresh the coordinates
   *  x[][] is the wrapped coordinates
   *  `n' may be less than x.length */
  public void refresh(double x[][], int n,
      boolean center, boolean adjScale) {
    if (model == null) {
      model = new XYZModel();
      adjScale = true;
    }
    model.updateXYZ(x, n, center);
    if ( adjScale ) realSpan = model.getSpan(x, n);
    repaint();
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
    tmpMatrix.xrot(360.0 * (mouseY - y) / getSize().height);
    tmpMatrix.yrot(360.0 * (x - mouseX) / getSize().width);
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


