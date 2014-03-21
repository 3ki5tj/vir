/** Computing the virial coefficients of the 2D hard sphere fluid
 *  by Monte Carlo simulation */
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



public class VirSampApp extends JApplet implements ActionListener
{
  MCSamp mc = new MCSamp(3, 6);
  int delay = 100; // interval of repainting
  int speed = 10000; // mc steps per frame
  Timer timer;

  XYScheme scheme; // schematic plot
  XYZCanvasVir canvas; // 3D animation

  JPanel cpnl, spnl, mpnl, mpnl2;
  JTextField      tDim     = new JTextField("" + mc.D);
  JTextField      tOrder   = new JTextField("" + mc.N);
  JTextField      tDelay   = new JTextField("" + delay);
  JTextField      tSpeed   = new JTextField("" + speed);
  JTextField      tMCAmp   = new JTextField("" + mc.mcAmp*mc.D);
  JTextField      tSampFreq = new JTextField("" + mc.sampFreq);
  JToggleButton   bStart   = new JToggleButton("Start", false);
  JButton         bReset   = new JButton("Reset");
  JButton         bClear   = new JButton("Clear data");
  JButton         bView    = new JButton("Reset view");
  JLabel          lStatus  = new JLabel("Status");
  JTextField      tVir     = new JTextField("0");
  JTextField      tAvStar  = new JTextField("0");
  JTextField      tAvRing  = new JTextField("0");
  JTextField      tAcc     = new JTextField("0");
  JCheckBox       chkMDS2D = new JCheckBox("Use 2D multi-dimensional scaling (MDS)");
  JCheckBox       chkMDS3D = new JCheckBox("Use 3D multi-dimensional scaling (MDS)");
  JCheckBox       chkMDS2DColor = new JCheckBox("Color according to 2D MDS");
  JCheckBox       chkMDS3DColor = new JCheckBox("Color according to 3D MDS");
  public static final long serialVersionUID = 1L;

  /** initialize the handler */
  public void init() {
    tDim.setHorizontalAlignment(JTextField.CENTER);
    tOrder.setHorizontalAlignment(JTextField.CENTER);
    tDelay.setHorizontalAlignment(JTextField.CENTER);
    tSpeed.setHorizontalAlignment(JTextField.CENTER);
    tMCAmp.setHorizontalAlignment(JTextField.CENTER);
    tSampFreq.setHorizontalAlignment(JTextField.CENTER);

    tVir.setHorizontalAlignment(JTextField.RIGHT);
    tAvStar.setHorizontalAlignment(JTextField.RIGHT);
    tAvRing.setHorizontalAlignment(JTextField.RIGHT);
    tAcc.setHorizontalAlignment(JTextField.RIGHT);

    Container box = getContentPane();
    box.setLayout(new BorderLayout());

    // cpnl is the panel for controls, such as buttons, textfields, etc.
    cpnl = new JPanel();
    cpnl.setLayout( new GridLayout(4, 6) ); // four rows, six columns
    box.add(cpnl, BorderLayout.NORTH);

    // add controls
    cpnl.add(bStart);
    bStart.addActionListener(this);

    cpnl.add(bReset);
    bReset.addActionListener(this);

    cpnl.add(new JLabel(" Dimension (D):"));
    tDim.addActionListener(this);
    cpnl.add(tDim);

    cpnl.add(new JLabel(" Order (n):"));
    tOrder.addActionListener(this);
    cpnl.add(tOrder);

    cpnl.add(new JLabel(" Delay (ms):"));
    tDelay.addActionListener(this);
    cpnl.add(tDelay);

    cpnl.add(new JLabel(" Steps/frame:"));
    tSpeed.addActionListener(this);
    cpnl.add(tSpeed);

    cpnl.add(new JLabel(" Samp. freq.:"));
    tSampFreq.addActionListener(this);
    cpnl.add(tSampFreq);

    cpnl.add(new JLabel(" Bn/B2^(n-1):"));
    tVir.setEditable(false);
    cpnl.add(tVir);

    cpnl.add(new JLabel(" < HSFb >:"));
    tAvStar.setEditable(false);
    cpnl.add(tAvStar);

    cpnl.add(new JLabel(" < Ring >:"));
    tAvRing.setEditable(false);
    cpnl.add(tAvRing);

    cpnl.add(new JLabel(" Move size D*\u0394X:"));
    tMCAmp.addActionListener(this);
    cpnl.add(tMCAmp);

    cpnl.add(new JLabel(" Acc. ratio:"));
    tAcc.setEditable(false);
    cpnl.add(tAcc);

    cpnl.add(bClear);
    bClear.addActionListener(this);

    cpnl.add(bView);
    bView.addActionListener(this);

    // create a panel for animation and schematic plot
    mpnl = new JPanel();
    mpnl.setLayout(new GridLayout(1, 2));
    box.add(mpnl, BorderLayout.CENTER);

    // scheme is the schematic diagram
    scheme = new XYScheme();
    scheme.setMC(mc);
    mpnl.add(scheme);

    // canvas is the place for 3D animation
    canvas = new XYZCanvasVir();
    mpnl.add(canvas);

    // a thin bar below `scheme' and `canvas'
    mpnl2 = new JPanel();
    mpnl2.setLayout(new GridLayout(2, 2));
    mpnl2.add(chkMDS2D);
    mpnl2.add(chkMDS3D);
    chkMDS2D.addActionListener(this);
    chkMDS3D.addActionListener(this);
    chkMDS2DColor.setEnabled(false);
    chkMDS3DColor.setEnabled(false);
    mpnl2.add(chkMDS2DColor);
    mpnl2.add(chkMDS3DColor);
    chkMDS2DColor.addActionListener(this);
    chkMDS3DColor.addActionListener(this);

    lStatus.setFont(new Font("Courier", Font.BOLD, 14));

    // spnl is the panel at the bottom
    spnl = new JPanel(); // create a panel for status
    spnl.setPreferredSize(new Dimension(800, 70));
    box.add(spnl, BorderLayout.SOUTH);

    // add mpnl2 and the status bar
    //spnl.setLayout(new GridLayout(2, 1));
    spnl.setLayout(new BorderLayout());
    spnl.add(mpnl2, BorderLayout.PAGE_START);
    spnl.add(lStatus, BorderLayout.PAGE_END);

    Font font = new Font("Arial", Font.BOLD, 14);
    tVir.setFont(font);

    timer = new Timer(delay, this);
    timer.start();
    timer.stop();
  }

  public void start() { }

  public void update(Graphics g) {
    paint(g);
  }

  /** Paint the applet */
  public void paint(Graphics g) {
    double vir = mc.getVir();
    lStatus.setText("steps = " + mc.step + ", "
         + "D = " + mc.D + ", "
         + "n = " + mc.N + ", "
         + "vir = " + String.format("%+.5e", vir) + ", "
         + "avFb = " + String.format("%+8.5f", mc.avFb.getAve()) + ", "
         + "avNr = " + String.format("%+8.5f", mc.avNr.getAve()) );
    tVir.setText(String.format("%.5e", vir) + "  ");
    tAvStar.setText(String.format("%.6f", mc.avFb.getAve()) + "  ");
    tAvRing.setText(String.format("%.6f", mc.avNr.getAve()) + "  ");
    tAcc.setText(String.format("%.4f", mc.mcTot > 0 ? 1.0*mc.mcAcc/mc.mcTot : 0.0) + "  ");
    scheme.repaint();
    mc.rmcom();
    canvas.refresh(mc.x, mc.N, mc.xo, mc.xi,
        false, mc.moveAtom, mc.moveAcc, false);
    cpnl.repaint();
    spnl.repaint();
  }

  /** Event handler */
  public void actionPerformed(ActionEvent e) {
    Object src = e.getSource();

    boolean adjScale = false;

    // the regular time for the regular animation
    if (src == timer) {
      for (int i = 0; i < speed; i++) // sample a few steps
        mc.mcStep();
      //mc.g.print(); // print the diagram
      repaint(); // calls the paint() function above
      return;
    }

    // the user may have changed the value of some parameters
    // but forget to press Enter to trigger an event
    // so we always refresh all parameters when an effect occurs
    int dim = mc.D;
    try {
      dim = Integer.parseInt(tDim.getText().trim());
    } catch (NumberFormatException err) {}
    if (dim < 2) {
      dim = 2;
      tDim.setText("" + dim);
    } else if (dim > 100) {
      dim = 100;
      tDim.setText("" + dim);
    }

    int n = mc.N;
    try {
      n = Integer.parseInt(tOrder.getText().trim());
    } catch (NumberFormatException err) {}
    if (n < 2) {
      n = 2;
      tOrder.setText("" + n);
    } else if (n > Diagram.NMAX) {
      n = Diagram.NMAX;
      tOrder.setText("" + n);
    }

    int delayOld = timer.getDelay();
    try {
      delay = Integer.parseInt(tDelay.getText().trim());
    } catch (NumberFormatException err) {}
    if (delay < 1) delay = 1;
    if (delay != delayOld)
      tDelay.setText("" + delay);

    try {
      speed = Integer.parseInt(tSpeed.getText().trim());
    } catch (NumberFormatException err) {}
    if (speed < 1) {
      speed = 1;
      tSpeed.setText("" + speed);
    }

    int sampFreq = mc.sampFreq;
    try {
      sampFreq = Integer.parseInt(tSampFreq.getText().trim());
    } catch (NumberFormatException err) {}
    mc.sampFreq = sampFreq;

    double amp = mc.mcAmp * mc.D;
    try {
      amp = Double.parseDouble(tMCAmp.getText().trim());
    } catch (NumberFormatException err) {}
    if (amp < 0.) {
      amp = -amp;
      tMCAmp.setText("" + amp);
    }
    mc.mcAmp = amp / dim;

    if (src == tDelay || delay != delayOld) {
      boolean running = timer.isRunning();
      if ( running ) timer.stop();
      timer = new Timer(delay, this);
      if ( running ) timer.start();
    }

    if (src == tSampFreq || src == tMCAmp
        || src == bReset || src == bClear) {
      mc.clearData();
    }

    if (src == bStart) {
      boolean on = bStart.isSelected();
      if (on) {
        timer.restart();
        //System.out.printf("delay %d\n", delay);
        bStart.setText("Pause");
      } else {
        timer.stop();
        bStart.setText("Resume");
      }
    }

    if (src == tDim || src == tOrder || src == bReset
     || dim != mc.D || n != mc.N) {
      if ( src == bReset ) {
        if ( timer.isRunning() ) timer.stop();
        bStart.setSelected(false);
        bStart.setText("Start");
      }
      mc.init(dim, n);
      mc.moveAtom = -1;
      adjScale = true;
    }

    if (src == bView) { // reset the view
      canvas.zoomScale = 1;
      canvas.viewMatrix.unit();
      adjScale = true;
    }

    // change the mode of multidimensional scaling
    if (src == chkMDS2D) {
      scheme.useMDS = chkMDS2D.isSelected();
      chkMDS2DColor.setEnabled( scheme.useMDS );
    }

    if (src == chkMDS2DColor && chkMDS2D.isSelected()) {
      scheme.colorMDS = chkMDS2DColor.isSelected();
    }

    scheme.repaint();

    if (src == chkMDS3D) {
      canvas.useMDS = chkMDS3D.isSelected();
      chkMDS3DColor.setEnabled( canvas.useMDS );
      adjScale = true;
    }

    if (src == chkMDS3DColor && chkMDS3D.isSelected()) {
      canvas.colorMDS = chkMDS3DColor.isSelected();
      adjScale = true;
    }

    mc.rmcom();
    canvas.refresh(mc.x, mc.N, mc.xo, mc.xi,
        false, mc.moveAtom, mc.moveAcc, adjScale);

    repaint();
  }
}



