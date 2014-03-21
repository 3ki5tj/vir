import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import static java.lang.Math.*;



/** Panel for drawing the schematic diagram */
class XYScheme extends JPanel
{
  Image img;
  Graphics imgG;
  Dimension imgSize;
  MCSamp mc;
  Atom atoms[];
  public static final long serialVersionUID = 2L;

  public XYScheme() {
    super();
    setBackground( Color.WHITE );
  }

  void setMC(MCSamp m) { mc = m; }

  public boolean newImgBuf() {
    Dimension sz = getSize();
    if (sz.width == 0 || sz.height == 0)
      return false;
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

  private static final Color colorLine         = Color.BLACK;
  private static final Color colorLargeCircle  = new Color(230, 230, 230);
  private static final Atom atomDefault = new Atom(0.1, 0.2, 1.0, 1.0, 0.6, 0.4, 1.0);
  private static final Atom atomFailed = new Atom(1.0, 0.1, 0.1, 1.0, 0.5, 0.3, 2.0);
  private static final Atom atomMoved  = new Atom(0.1, 1.0, 0.1, 1.0, 0.7, 0.3, 2.0);

  boolean useMDS = false;
  boolean colorMDS = false;

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    newImgBuf();
    Dimension sz = getSize();
    imgG.setColor( Color.WHITE );
    int w = getSize().width;
    int h = getSize().height;
    imgG.fillRect(0, 0, w, h);

    if (mc == null) return;
    int np = mc.N;

    // make lines thicker
    Graphics2D imgG2 = (Graphics2D) imgG;
    int lineWidth = (int) (40./np);
    if (lineWidth > 4) lineWidth = 4;
    else if (lineWidth < 1) lineWidth = 1;
    imgG2.setStroke(new BasicStroke(lineWidth));

    int rad = 400 / np; // radius of the balls
    if (rad > 40) rad = 40;
    int [][] xy;

    if (useMDS) {
      xy = getMDS(np, rad);
    } else {
      int radius = (int) (w/2 - rad - 5); // the big circle
      xy = new int [np][2];
      // compute the 2D coordinates around the circle
      for (int i = 0; i < np; i++) {
        double theta = 2*PI*i/np - PI*.5;
        xy[i][0] = (int) (w/2 + radius * cos(theta) + .5);
        xy[i][1] = (int) (h/2 + radius * sin(theta) + .5);
      }
      // draw the big circle
      imgG.setColor( colorLargeCircle );
      imgG.drawOval(w/2 - radius, h/2 - radius, 2*radius, 2*radius);
    }

    // draw connections
    imgG.setColor( colorLine );
    for (int i = 0; i < np; i++) {
      for (int j = i + 1; j < np; j++) {
        if ( mc.g.isLinked(i, j) ) {
          // draw a line between i and j
          imgG2.drawLine(xy[i][0], xy[i][1], xy[j][0], xy[j][1]);
        }
      }
    }

    // draw the small circles
    for (int i = 0; i < np; i++) {
      Atom atom = atomDefault;
      if (useMDS && colorMDS) atom = atoms[i];
      if (i == mc.moveAtom)
        atom = mc.moveAcc ? atomMoved : atomFailed;
      atom.paint(imgG, xy[i][0], xy[i][1], 15, rad);
    }

    g.drawImage(img, 0, 0, this);
  }

  /** Get the best coordinates from multidimensional scaling */
  int [][] getMDS(int np, int rad) {
    double [][] xy0 = new double [np][2];
    int [][] xy = new int [np][2];

    if (atoms == null || atoms.length != np)
      atoms = new Atom [np];

    for (int i = 0; i < np; i++) {
      double theta = 2*PI*i/np - PI*.5;
      xy0[i][0] = cos(theta);
      xy0[i][1] = sin(theta);
    }

    MDS mds = new MDS(mc.x);
    mds.min(xy0);
    if ( colorMDS ) mds.setAtoms(xy0, atoms);

    int w = getSize().width, h = getSize().height;
    int margin = 5 + rad;

    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    for (int i = 0; i < np; i++) {
      if (xy0[i][0] > xmax) xmax = xy0[i][0];
      else if (xy0[i][0] < xmin) xmin = xy0[i][0];
      if (xy0[i][1] > ymax) ymax = xy0[i][1];
      else if (xy0[i][1] < ymin) ymin = xy0[i][1];
    }
    xmax = max(xmax, -xmin);
    ymax = max(ymax, -ymin);

    double scl = min((.5 * w - margin) / xmax, (.5 * h - margin) / ymax);
    for (int i = 0; i < np; i++) {
      xy[i][0] = (int) (w / 2 + scl * xy0[i][0]);
      xy[i][1] = (int) (h / 2 + scl * xy0[i][1]);
    }
    xy0 = null;
    return xy;
  }
}



