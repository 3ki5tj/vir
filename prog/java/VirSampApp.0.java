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

  XYZCanvas canvas; // 3D animation
  XYScheme scheme; // schematic plot

  JPanel cpnl, spnl, mpnl, mpnl2;
  JTextField      tDim     = new JTextField("    " + mc.D);
  JTextField      tOrder   = new JTextField("    " + mc.N);
  JTextField      tDelay   = new JTextField("   " + delay);
  JTextField      tSpeed   = new JTextField("   " + speed);
  JTextField      tMCAmp   = new JTextField("   " + mc.mcAmp*mc.D);
  JTextField      tSampFreq = new JTextField("   " + mc.sampFreq);
  JToggleButton   bStart   = new JToggleButton("Start", false);
  JButton         bReset   = new JButton("Reset");
  JButton         bClear   = new JButton("Clear");
  JLabel          lStatus  = new JLabel("Status");
  JTextField      tVir     = new JTextField("   0");
  JTextField      tAvStar  = new JTextField("   0");
  JTextField      tAvRing  = new JTextField("   0");
  JTextField      tAcc     = new JTextField("   0");
  JCheckBox       chkMDS2D = new JCheckBox("Use multi-dimensional scaling (2D)");
  JCheckBox       chkMDS3D = new JCheckBox("Use multi-dimensional scaling (3D)");
  public static final long serialVersionUID = 1L;

  /** initialize the handler */
  public void init() {
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

    // create a panel for animation and schematic plot
    mpnl = new JPanel();
    mpnl.setLayout(new GridLayout(1, 2));
    box.add(mpnl, BorderLayout.CENTER);

    // scheme is the schematic diagram
    scheme = new XYScheme();
    scheme.setMC(mc);
    mpnl.add(scheme);

    // canvas is the place for 3D animation
    canvas = new XYZCanvas();
    mpnl.add(canvas);

    // a thin bar below `scheme' and `canvas'
    mpnl2 = new JPanel();
    mpnl2.setLayout(new GridLayout(1, 2));
    mpnl2.add(chkMDS2D);
    mpnl2.add(chkMDS3D);
    chkMDS2D.addActionListener(this);
    chkMDS3D.addActionListener(this);

    lStatus.setFont(new Font("Courier", Font.BOLD, 14));

    // spnl is the panel at the bottom
    spnl = new JPanel(); // create a panel for status
    box.add(spnl, BorderLayout.SOUTH);

    // add mpnl2 and the status bar
    spnl.setLayout(new GridLayout(2, 1));
    spnl.add(mpnl2, BorderLayout.NORTH);
    spnl.add(lStatus);

    Font font = new Font("Arial", Font.BOLD, 14);
    tVir.setFont(font);

    timer = new Timer(delay, this);
    timer.start();
    timer.stop();
  }

  public void start() {
    scheme.start();
    canvas.gauge(mc.x, mc.N);
    repaint();
  }

  public void update(Graphics g) {
    paint(g);
  }

  /** paint the applet */
  public void paint(Graphics g) {
    double vir = mc.getVir();
    lStatus.setText("steps = " + mc.step + ", "
         + "D = " + mc.D + ", "
         + "n = " + mc.N + ", "
         + "vir = " + String.format("%+.5e", vir) + ", "
         + "avfb = " + String.format("%+8.5f", mc.avfb.getAve()) + ", "
         + "avnr = " + String.format("%+8.5f", mc.avnr.getAve())
         + "");
    tVir.setText(" " + String.format("%.5e", vir));
    tAvStar.setText(" " + String.format("%g", mc.avfb.getAve()));
    tAvRing.setText(" " + String.format("%g", mc.avnr.getAve()));
    tAcc.setText(" " + String.format("%g", mc.mctot > 0 ? 1.0*mc.mcacc/mc.mctot : 0.0));
    scheme.repaint();
    if (canvas.model != null) {
      mc.rmcom();
      canvas.update(mc.x, mc.N,
                    mc.moveAcc ? mc.xo : mc.xi,
                    mc.moveAtom, mc.moveAcc);
      canvas.repaint();
    }
    cpnl.repaint();
    spnl.repaint();
  }

  /** Event handler */
  public void actionPerformed(ActionEvent e) {
    Object src = e.getSource();

    // the regular time for the regular animation
    if (src == timer) {
      for (int i = 0; i < speed; i++) // sample a few steps
        mc.mcStep();
      //mc.g.print();
      repaint(); // calls the paint() function above
      return;
    }

    //System.out.println( "trigger: " + src.toString() );

    // the user may have changed the value of some parameters
    // but forget to press Enter to trigger an event
    // so we always refresh all parameters when an effect occurs
    int dim = mc.D;
    try {
      dim = Integer.parseInt(tDim.getText().trim());
    } catch (NumberFormatException err) {}
    if (dim < 2) {
      dim = 2;
      tDim.setText(" " + dim);
    } else if (dim > 100) {
      dim = 100;
      tDim.setText(" " + dim);
    }

    int n = mc.N;
    try {
      n = Integer.parseInt(tOrder.getText().trim());
    } catch (NumberFormatException err) {}
    if (n < 2) {
      n = 2;
      tOrder.setText(" " + n);
    } else if (n > Diagram.NMAX) {
      n = Diagram.NMAX;
      tOrder.setText(" " + n);
    }

    int delay_old = timer.getDelay();
    try {
      delay = Integer.parseInt(tDelay.getText().trim());
    } catch (NumberFormatException err) {}
    if (delay < 1) delay = 1;
    if (delay != delay_old)
      tDelay.setText(" " + delay);

    try {
      speed = Integer.parseInt(tSpeed.getText().trim());
    } catch (NumberFormatException err) {}
    if (speed < 1) {
      speed = 1;
      tSpeed.setText(" " + speed);
    }

    int sampFreq = mc.sampFreq;
    try {
      sampFreq = Integer.parseInt(tSampFreq.getText().trim());
    } catch (NumberFormatException err) {}
    if (sampFreq < 1) {
      sampFreq = 1;
      tSpeed.setText(" " + sampFreq);
    }
    mc.sampFreq = sampFreq;

    double amp = mc.mcAmp * mc.D;
    try {
      amp = Double.parseDouble(tMCAmp.getText().trim());
    } catch (NumberFormatException err) {}
    if (amp < 0.) {
      amp = -amp;
      tMCAmp.setText(" " + amp);
    }
    mc.mcAmp = amp / dim;

    if (src == tDelay || delay != delay_old) {
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
        System.out.printf("delay %d\n", delay);
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
      canvas.viewMatrix.unit(); // reset the view
      canvas.bStarted = false;
    }

    // change the mode of multidimensional scaling
    if (src == chkMDS2D)
      scheme.useMDS = chkMDS2D.isSelected();

    if (src == chkMDS3D)
      canvas.useMDS = chkMDS3D.isSelected();

    // if the schematic diagram has not been initialized,
    // initialize it
    if ( !scheme.bStarted )
      scheme.newImgBuf();

    // if the canvas has not been initialized, initialize it
    if ( !canvas.bStarted )
      canvas.gauge(mc.x, mc.N);
    repaint();
  }
}




