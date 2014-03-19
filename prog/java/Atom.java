import java.applet.*;
import java.awt.*;
import java.awt.image.*;



/** Atom: a ball */
class Atom {
  private final static int R = 120;
  private final static int hx = 45; // (hx, hy) is the offset of the spot light from the center
  private final static int hy = 45;
  private static int maxr; // maximal distance from the spotlight
  public final static int nBalls = 16; // shades of grey
  // spotlight brightness (0.0, 1.0)
  // 1.0 means the spotlight is pure white
  // 0.0 means the spotlight is the same as the solid color
  double spotlightAmp = .4;
  // color damping along the distance from the spotlight
  double rContrast = 0.7;
  // z-depth contrast (0.0, inf)
  // inf means the fartherest atom is completely dark
  // 0.0 means the fartherest atom is the same as the
  //     nearest atom
  double zContrast = 2.0;

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
  private int Rl = 100, Gl = 100, Bl = 100;
  private Image balls[]; // 0..nBalls-1, at different z distances

  /** Constructor */
  Atom(int r, int g, int b) {
    setRGB(r, g, b);
  }

  /** colors that range from 0 to 1 */
  Atom(double r, double g, double b) {
    setRGB(r, g, b);
  }

  double relRadius = 1; // only used to save the radius

  Atom(double r, double g, double b, double rad) {
    this(r, g, b);
    relRadius = rad;
  }

  Atom(double r, double g, double b, double rad,
       double rcontrast, double spotlight, double zcontrast) {
    relRadius = rad;
    rContrast = rcontrast;
    spotlightAmp = spotlight;
    zContrast = zcontrast;
    setRGB(r, g, b);
  }

  /** Set color */
  void setRGB(int r, int g, int b) {
    Rl = r;
    Gl = g;
    Bl = b;
    makeBalls();
  }

  void setRGB(double r, double g, double b) {
    Rl = (int) (256*r - 1e-6);
    Gl = (int) (256*g - 1e-6);
    Bl = (int) (256*b - 1e-6);
    makeBalls();
  }

  /** Linearly interpolate colors */
  private int blend(double fg, double bg, double fgfactor) {
    return (int) (bg + (fg - bg) * fgfactor);
  }

  // need a component instance to call createImage()
  private static Component component = new Applet();

  /** Prepare ball images with different sizes */
  private void makeBalls() {
    balls = new Image[nBalls];
    byte red[]   = new byte[256];
    byte green[] = new byte[256];
    byte blue[]  = new byte[256];
    for (int id = 0; id < nBalls; id++) {
      // smaller `b' means closer to black
      // if id == 0 (fartherest from the viewer)
      //        b = 1/(1 + zContrast)
      //        the outer blend() gives a color close to bgGrey
      // if id == nBalls - 1 (closest to the viewer),
      //        b = 1, the outer blend() gives the color of
      //        the inner blend()
      double b = (zContrast*id/(nBalls - 1) + 1) / (zContrast + 1);
      for (int i = maxr; i >= 0; --i) {
        // closeness to the spotlight
        double q = 1 - 1. * i / maxr;
        // dampness of the color along the radius
        double p = 1 - rContrast * i / maxr;
        // contrast of the spotlight
        // if i == 0 (closest to the spotlight),
        //        d = 1.0 - spotlightAmp, the inner
        //        blend() gives a color close to 255
        //        (if spotlightAmp == 1).
        // if i == maxr (fartherest from the spotlight),
        //        d = 1.0, the inner blend() gives
        //        the foreground color, i.e., Rl, Gl, Bl
        // Thus, the inner blend() depends on the distance
        // from the spotlight, i == 0 means to be closest
        // to the spotlight
        double d = 1 - q * spotlightAmp;
        red[i]   = (byte) blend(blend(Rl*p, 255, d), 0, b);
        green[i] = (byte) blend(blend(Gl*p, 255, d), 0, b);
        blue[i]  = (byte) blend(blend(Bl*p, 255, d), 0, b);
      }
      // 256 color model
      IndexColorModel model = new IndexColorModel(
          8, maxr + 1, red, green, blue, 0);
      balls[id] = component.createImage(
          new MemoryImageSource(R * 2, R * 2, model, data, 0, R * 2) );
    }
  }

  /** Draw a ball at screen coordinate (x, y) with a ball index `id'
   *  (0, 0) represents the top-left corner
   *  `x', `y', `radius' are given in pixels
   *  the ball index (gray code) `id' can be 0 to 15 */
  void paint(Graphics gc, int x, int y, int id, double radius) {
    if (balls == null) makeBalls();
    Image img = balls[id]; // id = [0..15]

    int size = (int) (radius * 2 + .5);
    gc.drawImage(img, x - size/2, y - size/2, size, size, null);
    //System.out.println("" + x + " " + y + " " + id + " " + radius);
  }
}




