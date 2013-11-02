/* Computing the virial coefficients of the 2D hard sphere fluid
   by Monte Carlo simulation

   Copyright 2013 Cheng Zhang

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License can be found in
   <http://www.gnu.org/licenses/>.
*/
import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import static java.lang.Math.*;
import java.text.*;
import javax.swing.*;
import java.util.Hashtable;
import java.util.Random;



class VirSamp {
  public int D = 3; /* the dimension */
  public int N = 7; /* the order */
  Diagram g, ng;

  int step; // simulation steps
  double mcAmp = 1.5; // amplitude of randomly moving a particle
  long mcacc = 0, mctot = 0;
  double x[][], xi[];
  Ave star = new Ave(), ring = new Ave();
  double vir;

  Random rng = new Random();

  /** constructor */
  VirSamp(int dim, int order) {
    init(dim, order);
  }

  /** initialize system */
  public void init(int dim, int order) {
    int i, j;
    D = dim;
    N = order;
    mcAmp = 1.5/D;
    x = new double [N][D];
    for (i = 0; i < N; i++)
      for (j = 0; j < D; j++)
        x[i][j] = rng.nextGaussian() * 1.0/Math.sqrt(D*N);
    xi = new double [D];
    g = new Diagram(N, D, x);
    ng = new Diagram(N, D, x);
  }

  void clearData() {
    step = 0;
    mcacc = mctot = 0;
    star.clear();
    ring.clear();
  }

  /** a step of the Metropolis algorithm */
  public void mcstep() {
    int i, k;
    i = (int) (N * rng.nextDouble());
    /* displace particle i */
    for (k = 0; k < D; k++) {
      xi[k] = x[i][k];
      x[i][k] += (2 * rng.nextDouble() - 1) * mcAmp;
    }
    ng.fromCoordinates(D, x);
    mctot += 1;
    if ( ng.biconnected() ) {
      mcacc += 1;
    } else { /* recover */
      for (k = 0; k < D; k++)
        x[i][k] = xi[k];
    }
  }

  /** remove the center of mass motion */
  void rmcom() {
    for (int k = 0; k < D; k++) {
      xi[k] = 0;
      for (int i = 0; i < N; i++)
        xi[k] += x[i][k];
      xi[k] /= N;
      for (int i = 0; i < N; i++)
        x[i][k] -= xi[k];
    }
  }

  /** get three dimensional coordinates */
  public double [] getxyz() {
    double [] xyz;
    int i, j;

    rmcom();
    xyz = new double [N*3];
    for (i = 0; i < N; i++) {
      xyz[i*3 + 0] = x[i][0];
      xyz[i*3 + 1] = x[i][1];
      if (D == 2) {
        xyz[i*3 + 2] = 0;
      } else {
        xyz[i*3 + 2] = x[i][2];
      }
    }
    return xyz;
  }

  private void rvn_diff(double z[], double x[], double y[]) {
    for (int i = 0; i < D; i++)
      z[i] = x[i] - y[i];
  }
}



class Diagram {
  public static final int nvmax = 62; /* maximal number of vertices */
  int n;
  long c[];

  /** construct a graph */
  Diagram(int nv) {
    if (nv > nvmax) nv = nvmax;
    n = nv;
    c = new long[n];
  }

  /** construct from D-dimensional coordinates */
  Diagram(int nv, int dim, double r[][]) {
    this(nv);
    fromCoordinates(dim, r);
  }

  /** update the graph form D-dimensional coordinates */
  void fromCoordinates(int dim, double r[][]) {
    int i, j, k;
    double r2, dx;
    empty();
    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++) {
        r2 = 0;
        for (k = 0; k < dim; k++) {
          dx = r[i][k] - r[j][k];
          r2 += dx * dx;
        }
        if (r2 < 1) link(i, j);
      }
    }
  }

  /** link two vertices i and j */
  void link(int i, int j) {
    c[i] |= 1L << j;
    c[j] |= 1L << i;
  }

  /** unlink two vertices i and j */
  void unlink(int i, int j) {
    c[i] &= ~(1L << j);
    c[j] &= ~(1L << i);
  }

  /** remove all edges in the graph */
  void empty() {
    for (int i = 0; i < n; i ++)
      c[i] = 0;
  }

  /** add all edges to the graph */
  void full() {
    long mask = mkbitsmask(n);
    for (int i = 0; i < n; i++) {
      c[i] = mask ^ mkbit(i);
    }
  }

  /** check if the subgraph of `vs' is connected */
  boolean connectedvs(long vs) {
    long stack = vs & (-vs), bk;
    int k;
    while (stack != 0) {
      bk = stack & (-stack);
      k = bit2id(bk);
      vs ^= bk; /* remove k from the to-do list */
      stack = (stack | c[k]) & vs;
      if (stack == vs)
        return true;
    }
    return false;
  }

  /** check if the graph is connected */
  boolean connected() {
    return connectedvs(mkbitsmask(n));
  }

  /** check if the graph is biconnected */
  boolean biconnected() {
    if (n > 2) {
      long mask = mkbitsmask(n);
      for (long b = 1; (b & mask) != 0; b <<= 1)
        if ( !connectedvs(mask ^ b) )
          return false;
      return true;
    } else if (n == 2) {
      return ((c[1] & 0x1l) != 0) ? true : false;
    } else return true;
  }

  /* below are bit operations */

  /** return a long integer with only ith lowest bit being 1 */
  long mkbit(int i) {
    return (long) 1 << i;
  }

  /** return a long integer with the lowest i bits being 1 */
  long mkbitsmask(int i) {
    return ((long) 1 << i) - 1;
  }

  private static final int [] bitrem = {
    -1,  0,  1, 39,  2, 15, 40, 23,  3, 12, 16, 59, 41, 19, 24, 54,
     4, -1, 13, 10, 17, 62, 60, 28, 42, 30, 20, 51, 25, 44, 55, 47,
     5, 32, -1, 38, 14, 22, 11, 58, 18, 53, 63,  9, 61, 27, 29, 50,
    43, 46, 31, 37, 21, 57, 52,  8, 26, 49, 45, 36, 56,  7, 48, 35,
     6, 34, 33};

  int bit2id(long bit) {
    return bitrem[(int) (bit % 67)];
  }

  int bitfirst(long x) {
    long bit = x & (-x);
    return bit2id(bit);
  }

  private static final int [] bytecount = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

  int bitcount(long x) {
    return bytecount[(int) (x         & 0xffl)] +
           bytecount[(int) ((x >>  8) & 0xffl)] +
           bytecount[(int) ((x >> 16) & 0xffl)] +
           bytecount[(int) ((x >> 24) & 0xffl)] +
           bytecount[(int) ((x >> 32) & 0xffl)] +
           bytecount[(int) ((x >> 40) & 0xffl)] +
           bytecount[(int) ((x >> 48) & 0xffl)];
  }
}



/* average and standard deviation */
class Ave {
  double cnt, xsm, x2sm;

  void clear() {
    cnt = xsm = x2sm = 0;
  }

  void add(double x) {
    cnt += 1;
    xsm += x;
    x2sm += x*x;
  }

  double getAve() {
    return (cnt > 0) ? xsm/cnt : 0;
  }

  double getStd() {
    if (cnt <= 1) return 0;
    double av = xsm/cnt;
    return sqrt(x2sm/cnt - av*av);
  }
}




public class VirSampApp extends JApplet implements ActionListener
{
  VirSamp mc = new VirSamp(3, 7);
  int delay = 100;
  int speed = 1000; // mc steps per frame
  Timer timer;

  MyCanvas canvas; // 3D drawing here

  JPanel cpnl, spnl;
  JTextField      tDim     = new JTextField("    " + mc.D);
  JTextField      tOrder   = new JTextField("    " + mc.N);
  JTextField      tSpeed   = new JTextField("   " + speed);
  JTextField      tMCAmp   = new JTextField("   " + mc.mcAmp);
  JToggleButton   bStart   = new JToggleButton("Start", false);
  JButton         bReset   = new JButton("Reset");
  JButton         bRetime  = new JButton("Reset time");
  JLabel          lStatus  = new JLabel("Status");
  JTextField      tVir     = new JTextField("   0");
  JTextField      tAvStar  = new JTextField("   0");
  JTextField      tAvRing  = new JTextField("   0");
  JTextField      tAcc     = new JTextField("   0");
  public static final long serialVersionUID = 1L;

  public void init() {
    Container box = getContentPane();
    box.setLayout(new BorderLayout());

    cpnl = new JPanel(); // create a panel for controls
    cpnl.setLayout(new GridLayout(19, 2));
    box.add(cpnl, BorderLayout.EAST);

    // add controls
    cpnl.add(bStart);
    bStart.addActionListener(this);

    cpnl.add(bReset);
    bReset.addActionListener(this);

    cpnl.add(new JLabel(" D:"));
    tDim.addActionListener(this);
    cpnl.add(tDim);

    cpnl.add(new JLabel(" N:"));
    tOrder.addActionListener(this);
    cpnl.add(tOrder);

    cpnl.add(new JLabel(" Steps/frame:"));
    tSpeed.addActionListener(this);
    cpnl.add(tSpeed);

    cpnl.add(new JLabel(" Move \u0394X:"));
    tMCAmp.addActionListener(this);
    cpnl.add(tMCAmp);

    cpnl.add(new JLabel(" Bn/B2^(n-1):"));
    tVir.setEditable(false);
    cpnl.add(tVir);

    cpnl.add(new JLabel(" < Star >:"));
    tAvStar.setEditable(false);
    cpnl.add(tAvStar);

    cpnl.add(new JLabel(" < Ring >:"));
    tAvRing.setEditable(false);
    cpnl.add(tAvRing);

    cpnl.add(new JLabel(" Acc. ratio:"));
    tAcc.setEditable(false);
    cpnl.add(tAcc);

    cpnl.add(bRetime);
    bRetime.addActionListener(this);

    spnl = new JPanel(); // create a panel for status
    box.add(spnl, BorderLayout.SOUTH);
    lStatus.setFont(new Font("Courier", 0, 12));
    spnl.add(lStatus);

    canvas = new MyCanvas();
    canvas.setmc(mc);
    canvas.preferredWidth = 2.0;
    box.add(canvas, BorderLayout.CENTER);

    timer = new Timer(delay, this);
    timer.start(); timer.stop();
  }

  public void start() {
    Atom.setApplet(this);
    canvas.start();
  }

  public void update(Graphics g) { paint(g); }

  DecimalFormat df = new DecimalFormat("###0.000");
  //int frame = 0;
  public void paint(Graphics g) {
    //System.out.println("frame: " + (frame++));
    lStatus.setText("n = " + mc.step + ", "
         + "N = " + mc.N + ", "
         + "vir = " + df.format(mc.vir) + ", "
         + "star = " + df.format(mc.star.getAve()) + ", "
         + "ring = " + df.format(mc.ring.getAve())
         + ";");
    tVir.setText(" " + df.format(mc.vir));
    tAvStar.setText(" " + df.format(mc.star.getAve()));
    tAvRing.setText(" " + df.format(mc.ring.getAve()));
    tAcc.setText(" " + df.format(mc.mctot > 0 ? 1.0*mc.mcacc/mc.mctot : 0.0));
    if (canvas.mdl != null) {
      canvas.mdl.updatexyz( mc.getxyz() );
      canvas.repaint();
    } else return;
    cpnl.repaint();
    spnl.repaint();
  }

  public void actionPerformed(ActionEvent e) {
    Object src = e.getSource();
    if (src == timer) {
      for (int i = 0; i < speed; i++) // sample a few steps
        mc.mcstep();
      repaint();
      return;
    }

    if (src == tDim || src == tOrder || src == bReset) {
      int dim = Integer.parseInt(tDim.getText().trim());
      int order = Integer.parseInt(tOrder.getText().trim());
      mc.init(dim, order);
      mc.clearData();
      canvas.prepColor();
    }

    if (src == tSpeed || src == bReset) {
      speed = Integer.parseInt(tSpeed.getText().trim());
      if (speed < 1) { speed = 1; tSpeed.setText("   " +speed); }
    }

    if (src == tMCAmp || src == bReset) {
      double amp = Double.parseDouble(tMCAmp.getText().trim());
      if (amp < 0.) { amp = -amp; tMCAmp.setText("   " + amp); }
      mc.mcAmp = amp;
      mc.clearData();
    }

    if (src == bRetime) {
      mc.clearData();
    }

    if (src == bStart) {
      boolean on = bStart.isSelected();
      if (on) {
        timer.restart();
        bStart.setText("Pause");
      } else {
        timer.stop();
        bStart.setText("Resume");
      }
    }

    if (src == tOrder || src == tDim) {
      int dim = Integer.parseInt(tDim.getText().trim());
      if (dim < 2) { dim = 2; tDim.setText(" " + dim); }
      int n = Integer.parseInt(tOrder.getText().trim());
      if (n < 2) { n = 2; tOrder.setText(" " + n); }
      mc.D = dim;
      mc.N = n;
      mc.init(mc.D, mc.N);
      canvas.preferredWidth = Math.pow(mc.N, 1./mc.D);
      canvas.gauge(); //start();
    }

    if (src == bReset) {
      if (timer.isRunning()) timer.stop();
      bStart.setSelected(false);
      bStart.setText("Start");
      mc.init(mc.D, mc.N);
    }
    repaint();
  }
}



/** panel for drawing particals */
class MyCanvas extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener {
  Image img;
  Graphics imgG;
  Dimension imgSize;

  double xfac;
  double xtheta, ytheta;
  double zoomscale = 1.0f;
  Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
  int prevx, prevy;
  String message;
  XYZModel mdl;
  VirSamp mc;
  public static final long serialVersionUID = 2L;

  public MyCanvas() {
    super();
    amat.yrot(2); amat.xrot(2);
  }

  void setmc(VirSamp m) { mc = m; }

  double preferredWidth = 0;

  /** determine scaling factors based on the given width */
  void gaugeByWidth(double xw) {
    if (mdl == null)
      mdl = new XYZModel(mc.getxyz());
    mdl.xmin = -xw * 0.5; mdl.xmax = xw * 0.5;
    mdl.ymin = -xw * 0.5; mdl.ymax = xw * 0.5;
    mdl.zmin = -xw * 0.5; mdl.zmax = xw * 0.5;

    double f1 = getSize().width / xw;
    double f2 = getSize().height / xw;
    xfac = 0.7f * (f1 < f2 ? f1 : f2);
    //System.out.println("xfac " + xfac + " xw " + xw + " preferred " + preferredWidth);
    newImgBuf();
    prepColor();
  }

  /** determine scaling factors based on the current coordinates */
  void gauge() {
    if (preferredWidth > 0) {
      gaugeByWidth(preferredWidth);
      return;
    }

    mdl = new XYZModel(mc.getxyz());
    /* find the largest span in three dimensions */
    mdl.findbox();
    double xw = mdl.xmax - mdl.xmin;
    double yw = mdl.ymax - mdl.ymin;
    double zw = mdl.zmax - mdl.zmin;
    if (yw > xw) xw = yw;
    if (zw > xw) xw = zw;
    System.out.println("canvas init: " + getSize().toString() + " " + xw + " " + yw + " " +zw);
    gaugeByWidth(xw);
  }

  public void start() {
    gauge();
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
  }

  /** set color of the default atom */
  void prepColor() {
    Atom.setRGB(0, 0, 255);
  }

  private void newImgBuf() {
    if (getSize().width == 0 || getSize().height == 0) {
      System.out.println("canvas not ready.\n");
      return;
    }
    //System.out.println("canvas " + getSize().toString());
    img = createImage(getSize().width, getSize().height);
    if (imgG != null) { imgG.dispose(); }
    imgG = img.getGraphics();
    imgSize = getSize();
  }

  public void update(Graphics g) {
    if (img == null)
      g.clearRect(0, 0, getSize().width, getSize().height);
    paintComponent(g);
  }

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    if (mdl != null) {
      mdl.mat.unit();
      mdl.mat.translate(-(mdl.xmin + mdl.xmax) / 2,
                        -(mdl.ymin + mdl.ymax) / 2,
                        -(mdl.zmin + mdl.zmax) / 2);
      mdl.mat.mult(amat);
      /* xfac is the default scaling from real coordinate
       * to the screen coordinats, it must be multiplied by
       * the zoom scale to get the actual scaling */
      double xfac1 = xfac * zoomscale;
      mdl.mat.scale(xfac1, xfac1, 16 * xfac1 / getSize().width);
      mdl.mat.translate(getSize().width / 2, getSize().height / 2, 8);
      mdl.xfac = xfac1; /* tell the XYZModel the scaling factor */
      //System.out.println("xfac " + xfac + " xfac1 " + xfac1 + " zoom " + zoomscale);
      mdl.transformed = false;
      if (img != null) {
        if (!imgSize.equals(getSize()))
          newImgBuf();
        imgG.setColor(Color.BLACK);
        imgG.fillRect(0, 0, getSize().width, getSize().height);
        mdl.paint(imgG);
        g.drawImage(img, 0, 0, this);
      } else
        mdl.paint(g);
    } else if (message != null) {
      g.drawString("Error in model:", 3, 20);
      g.drawString(message, 10, 40);
    }
  }

  /* event handling */
  public void mouseClicked(MouseEvent e) { }
  public void mousePressed(MouseEvent e) {
    prevx = e.getX();
    prevy = e.getY();
    e.consume();
  }
  public void mouseReleased(MouseEvent e) { }
  public void mouseEntered(MouseEvent e) { }
  public void mouseExited(MouseEvent e) { }

  public void mouseDragged(MouseEvent e) {
    int x = e.getX();
    int y = e.getY();
    tmat.unit();
    double xtheta = (prevy - y) * (360.0f / getSize().height);
    double ytheta = (x - prevx) * (360.0f / getSize().width);
    tmat.xrot(xtheta);
    tmat.yrot(ytheta);
    amat.mult(tmat);
    repaint();
    prevx = x;
    prevy = y;
    e.consume();
  }

  public void mouseMoved(MouseEvent e) { }

  public void mouseWheelMoved(MouseWheelEvent e) {
    int notches = e.getWheelRotation();
    if ((zoomscale -= 0.05f*notches) < 0.09999f) zoomscale = 0.1f;
    repaint();
  }
}


/** a set of atoms */
class XYZModel {
  double realxyz[]; /* 3D real coordinates 3*npmax */
  int screenxyz[];  /* 3D screen coordinates in pixels
                       only the first two dimensions are used in drawing */
  int zorder[]; /* z-order */
  int np, npmax;

  boolean transformed;
  Matrix3D mat; /* rotation/scaling/translation matrix for the conversion
                   from realxyz[] to screenxyz[] */
  double xfac = 1; /* real size --> screen size, to be set by MyCanvas */
  double ballSize = 0.5f; /* ball size (radius) in terms of the real coordinates
                             it is 0.5 because we are simulating hard spheres */
  double xmin, xmax, ymin, ymax, zmin, zmax; /* box dimension */

  /** constructor from the real coordinates */
  XYZModel(double r[]) {
    int n = r.length/3;
    mat = new Matrix3D();
    for (int i = 0; i < n; i++)
      addPoint(r[3*i], r[3*i+1], r[3*i+2]);
  }

  /** Add a point to this model */
  int addPoint(double x, double y, double z) {
    int i = np;
    if (i >= npmax)
      if (realxyz == null) {
        npmax = 100;
        realxyz = new double[npmax * 3];
      } else {
        npmax *= 2;
        double newxyz[] = new double[npmax * 3];
        System.arraycopy(realxyz, 0, newxyz, 0, realxyz.length);
        realxyz = newxyz;
      }
    i *= 3;
    realxyz[i]     = x;
    realxyz[i + 1] = y;
    realxyz[i + 2] = z;
    return np++;
  }

  /** Refresh coordinates */
  void updatexyz(double r[]) {
    int n = r.length/3;
    if (n > np) realxyz = new double[n * 3];
    np = n;
    for (int i = 0; i < 3*np; i++)
      realxyz[i] = r[i];
  }

  /** Transform all the points in this model */
  void transform() {
    if (transformed || np <= 0) return;
    if (screenxyz == null || screenxyz.length < np * 3)
      screenxyz = new int[np * 3];
    mat.transform(realxyz, screenxyz, np);
    transformed = true;
  }

  /** Paint this model to a graphics context.  It uses the matrix associated
   with this model to map from model space to screen space.
   The next version of the browser should have double buffering,
   which will make this *much* nicer */
  void paint(Graphics g) {
    if (realxyz == null || np <= 0) return;
    transform();
    if (zorder == null || zorder.length < np) {
      zorder = new int[np];
      for (int i = 0; i < np; i++) zorder[i] = i;
    }

    // bubble sort z-order
    for (int i = np - 1; --i >= 0;) {
      boolean flipped = false;
      for (int j = 0; j <= i; j++) {
        int a = zorder[j];
        int b = zorder[j+1];
        if (screenxyz[3*a + 2] > screenxyz[3*b + 2]) {
          zorder[j + 1] = a;
          zorder[j] = b;
          flipped = true;
        }
      }
      if (!flipped) break; // nothing after a sweep
    }

    double radius = xfac * ballSize;
    //System.out.println("xfac " + xfac + ", ball " + ballSize);
    for (int i = 0; i < np; i++) {
      int j = zorder[i];
      int grey = screenxyz[3*j + 2]; // adjust color according to z
      if (grey < 3) grey = 3;
      else if (grey > 15) grey = 15;
      /* screenxyz[3*j] and screenxyz[3*j + 1] are the (x, y) coordinates on screen */
      Atom.paint(g, screenxyz[3*j], screenxyz[3*j + 1], grey, radius);
    }
  }

  /** Find the bounding box of this model */
  void findbox() {
    if (np <= 0) return;
    double xmin = realxyz[0], xmax = xmin;
    double ymin = realxyz[1], ymax = ymin;
    double zmin = realxyz[2], zmax = zmin;
    for (int i = np * 3; (i -= 3) > 0;) {
      double x = realxyz[i], y = realxyz[i+1], z = realxyz[i+2];
      if (x < xmin) xmin = x; else if (x > xmax) xmax = x;
      if (y < ymin) ymin = y; else if (y > ymax) ymax = y;
      if (z < zmin) zmin = z; else if (z > zmax) zmax = z;
    }
    this.xmax = xmax; this.xmin = xmin;
    this.ymax = ymax; this.ymin = ymin;
    this.zmax = zmax; this.zmin = zmin;
  }
}



class Atom {
  private static Component applet;
  private static byte[] data;
  private final static int R = 40;
  private final static int hx = 15;  // (hx, hy) is the offset of the spot light from the center
  private final static int hy = 15;
  private final static int bgGrey = 0; // further atoms are to be blended with this color
  private final static int nBalls = 16;
  private final static double spotlightmag = .5f; // spotlight brightness, 1.f for full
  private static int maxr; // maximal intensity

  private static int Rl = 30, Gl = 10, Bl = 255;
  private static Image balls[]; // 0..nBalls-1, at different z distances
  public double ballSize; // ball size

  static { /* data[] is a bitmap image of the ball of radius R */
    data = new byte[R * 2 * R * 2];
    int mr = 0;
    for (int Y = 2 * R; --Y >= 0;) {
      int x0 = (int) (Math.sqrt(R * R - (Y - R) * (Y - R)) + 0.5);
      int p = Y * (R * 2) + R - x0;
      for (int X = -x0; X < x0; X++) {
        int x = X + hx;
        int y = Y - R + hy;
        int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
        if (r > mr) mr = r;
        data[p++] = r <= 0 ? 1 : (byte) r;
      }
    }
    maxr = mr;
  }

  /** set color */
  static void setRGB(int r, int g, int b) {
    Rl = r; Gl = g; Bl = b;
    mkballs();
  }

  static void setApplet(Component app) { applet = app; }

  private final static int blend(int fg, int bg, double fgfactor)
    { return (int) (bg + (fg - bg) * fgfactor); }

  /** prepare ball images with different sizes */
  private static void mkballs() {
    balls = new Image[nBalls];
    byte red[] = new byte[256];
    red[0] = (byte) bgGrey;
    byte green[] = new byte[256];
    green[0] = (byte) bgGrey;
    byte blue[] = new byte[256];
    blue[0] = (byte) bgGrey;
    for (int r = 0; r < nBalls; r++) {
      double b = (r + 1.) / nBalls;
      for (int i = maxr; i >= 1; --i) {
        double d = 1. * i / maxr;
        d = 1 - (1 - d) * spotlightmag;
        red[i] = (byte) blend(blend(Rl, 255, d), bgGrey, b);
        green[i] = (byte) blend(blend(Gl, 255, d), bgGrey, b);
        blue[i] = (byte) blend(blend(Bl, 255, d), bgGrey, b);
      }
      // 256 color model
      IndexColorModel model = new IndexColorModel(8, maxr + 1, red, green, blue, 0);
      balls[r] = applet.createImage(new MemoryImageSource(R * 2, R * 2, model, data, 0, R * 2));
    }
  }

  /** Draw a ball at screen coordinate (x, y) with a ball index `r' and
   *  the radius being `radius'. (0, 0) represents the top-left corner
   *  x, y, radius are given in terms of pixels
   *  the ball index (gray code) can be 0 to 15 */
  static void paint(Graphics gc, int x, int y, int r, double radius) {
    if (balls == null) mkballs();
    Image img = balls[r]; // r = [0..15]

    radius *= (36.0 + r) / (36 + 15); // make remote balls smaller
    int size = (int) (radius*2);
    gc.drawImage(img, x - size/2, y - size/2, size, size, null);
    //System.out.println("" + x + " " + y + " " + r + " " + radius);
  }
}



class Matrix3D {
  double xx, xy, xz, xo;
  double yx, yy, yz, yo;
  double zx, zy, zz, zo;
  static final double pi = 3.14159265;

  /** Create a new unit matrix */
  Matrix3D() {
    xx = 1.0f;
    yy = 1.0f;
    zz = 1.0f;
  }

  /** Scale along each axis independently */
  void scale(double xf, double yf, double zf) {
    xx *= xf; xy *= xf; xz *= xf; xo *= xf;
    yx *= yf; yy *= yf; yz *= yf; yo *= yf;
    zx *= zf; zy *= zf; zz *= zf; zo *= zf;
  }

  /** Translate the origin */
  void translate(double x, double y, double z) { xo += x; yo += y; zo += z; }

  /** rotate theta degrees about the y axis */
  void yrot(double theta) {
    theta *= (pi / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nxx = (xx * ct + zx * st);
    double Nxy = (xy * ct + zy * st);
    double Nxz = (xz * ct + zz * st);
    double Nxo = (xo * ct + zo * st);

    double Nzx = (zx * ct - xx * st);
    double Nzy = (zy * ct - xy * st);
    double Nzz = (zz * ct - xz * st);
    double Nzo = (zo * ct - xo * st);

    xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
    zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
  }

  /** rotate theta degrees about the x axis */
  void xrot(double theta) {
    theta *= (pi / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (yx * ct + zx * st);
    double Nyy = (yy * ct + zy * st);
    double Nyz = (yz * ct + zz * st);
    double Nyo = (yo * ct + zo * st);

    double Nzx = (zx * ct - yx * st);
    double Nzy = (zy * ct - yy * st);
    double Nzz = (zz * ct - yz * st);
    double Nzo = (zo * ct - yo * st);

    yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
    zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
  }

  /** rotate theta degrees about the z axis */
  void zrot(double theta) {
    theta *= pi / 180;
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (yx * ct + xx * st);
    double Nyy = (yy * ct + xy * st);
    double Nyz = (yz * ct + xz * st);
    double Nyo = (yo * ct + xo * st);

    double Nxx = (xx * ct - yx * st);
    double Nxy = (xy * ct - yy * st);
    double Nxz = (xz * ct - yz * st);
    double Nxo = (xo * ct - yo * st);

    yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
    xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
  }

  /** Multiply this matrix by a second: M = M*R */
  void mult(Matrix3D rhs) {
    double lxx = xx * rhs.xx + yx * rhs.xy + zx * rhs.xz;
    double lxy = xy * rhs.xx + yy * rhs.xy + zy * rhs.xz;
    double lxz = xz * rhs.xx + yz * rhs.xy + zz * rhs.xz;
    double lxo = xo * rhs.xx + yo * rhs.xy + zo * rhs.xz + rhs.xo;

    double lyx = xx * rhs.yx + yx * rhs.yy + zx * rhs.yz;
    double lyy = xy * rhs.yx + yy * rhs.yy + zy * rhs.yz;
    double lyz = xz * rhs.yx + yz * rhs.yy + zz * rhs.yz;
    double lyo = xo * rhs.yx + yo * rhs.yy + zo * rhs.yz + rhs.yo;

    double lzx = xx * rhs.zx + yx * rhs.zy + zx * rhs.zz;
    double lzy = xy * rhs.zx + yy * rhs.zy + zy * rhs.zz;
    double lzz = xz * rhs.zx + yz * rhs.zy + zz * rhs.zz;
    double lzo = xo * rhs.zx + yo * rhs.zy + zo * rhs.zz + rhs.zo;

    xx = lxx; xy = lxy; xz = lxz; xo = lxo;
    yx = lyx; yy = lyy; yz = lyz; yo = lyo;
    zx = lzx; zy = lzy; zz = lzz; zo = lzo;
  }

  /** Reinitialize to the unit matrix */
  void unit() {
    xo = 0; xx = 1; xy = 0; xz = 0;
    yo = 0; yx = 0; yy = 1; yz = 0;
    zo = 0; zx = 0; zy = 0; zz = 1;
  }

  /** Transform np points from v into tv.  v contains the input
   coordinates in floating point.  Three successive entries in
   the array constitute a point.  tv ends up holding the transformed
   points as integers; three successive entries per point */
  void transform(double v[], int tv[], int np) {
    double lxx = xx, lxy = xy, lxz = xz, lxo = xo;
    double lyx = yx, lyy = yy, lyz = yz, lyo = yo;
    double lzx = zx, lzy = zy, lzz = zz, lzo = zo;
    for (int i = np * 3; (i -= 3) >= 0;) {
      double x = v[i], y = v[i + 1], z = v[i + 2];
      tv[i]     = (int) (x * lxx + y * lxy + z * lxz + lxo);
      tv[i + 1] = (int) (x * lyx + y * lyy + z * lyz + lyo);
      tv[i + 2] = (int) (x * lzx + y * lzy + z * lzz + lzo);
    }
  }
}
