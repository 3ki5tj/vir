import java.applet.*;
import java.awt.*;
import java.awt.image.*;



/** Atom: a ball */
class Atom {
  private final static int R = 120;
  private final static int hx = 45;  // (hx, hy) is the offset of the spot light from the center
  private final static int hy = 45;
  private final static int bgGrey = 0; // further atoms are to be blended with this color
  private final static int nBalls = 16;
  private final static double spotlightmag = .9f; // spotlight brightness, 1.f for full
  private static int maxr; // maximal intensity

  private static byte data[];
  static { // data[] is a bitmap image of the ball of radius R
    data = new byte[R * 2 * R * 2];
    for (int Y = -R; Y < R; Y++) {
      int x0 = (int) (Math.sqrt(R * R - Y * Y) + 0.5);
      for (int X = -x0; X < x0; X++) {
        // sqrt(x^2 + y^2) gives distance from the spot light
        int x = X + hx, y = Y + hy;
        int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
        // set the maximal intensity to the maximal distance
        // (in pixels) from the spot light
        if (r > maxr) maxr = r;
        data[(Y + R) * (R * 2) + (X + R)] = (r <= 0) ? 1 : (byte) r;
      }
    }
  }

  // the following variables are atom dependent
  private int Rl = 30, Gl = 10, Bl = 255;
  private Image balls[]; // 0..nBalls-1, at different z distances
  public double ballSize; // ball size

  /** Constructor */
  Atom(int r, int g, int b) {
    setRGB(r, g, b);
  }

  /** Set color */
  void setRGB(int r, int g, int b) {
    Rl = r; Gl = g; Bl = b;
    makeBalls();
  }

  /** Linearly interpolate colors */
  private int blend(int fg, int bg, double fgfactor) {
    return (int) (bg + (fg - bg) * fgfactor);
  }

  // need a component instance to call createImage()
  private static Component component = new Applet();

  /** Prepare ball images with different sizes */
  private void makeBalls() {
    balls = new Image[nBalls];
    byte red[] = new byte[256];
    red[0] = (byte) bgGrey;
    byte green[] = new byte[256];
    green[0] = (byte) bgGrey;
    byte blue[] = new byte[256];
    blue[0] = (byte) bgGrey;
    for (int r = 0; r < nBalls; r++) {
      // smaller b means greyer
      double b = (r + 1. + 7) / (nBalls + 7);
      for (int i = maxr; i >= 1; --i) {
        // contrast of the spotlight
        double d = 1 - (1 - 1. * i / maxr) * spotlightmag;
        red[i]    = (byte) blend(blend(Rl, 255, d), bgGrey, b);
        green[i]  = (byte) blend(blend(Gl, 255, d), bgGrey, b);
        blue[i]   = (byte) blend(blend(Bl, 255, d), bgGrey, b);
      }
      // 256 color model
      IndexColorModel model = new IndexColorModel(
          8, maxr + 1, red, green, blue, 0);
      balls[r] = component.createImage(
          new MemoryImageSource(R * 2, R * 2, model, data, 0, R * 2) );
    }
  }

  /** Draw a ball at screen coordinate (x, y) with a ball index `r'
   *  the `radius' is measured in pixels
   *  (0, 0) represents the top-left corner
   *  x, y, radius are given in terms of pixels
   *  the ball index (gray code) can be 0 to 15 */
  void paint(Graphics gc, int x, int y, int r, double radius) {
    if (balls == null) makeBalls();
    Image img = balls[r]; // r = [0..15]

    int size = (int) (radius * 2 + .5);
    gc.drawImage(img, x - size/2, y - size/2, size, size, null);
    //System.out.println("" + x + " " + y + " " + r + " " + radius);
  }
}




