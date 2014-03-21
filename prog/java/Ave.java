import java.lang.Math.*;



/** Average and standard deviation */
class Ave {
  double cnt, xsm, x2sm;

  public void clear() {
    cnt = xsm = x2sm = 0;
  }

  public void add(double x) {
    cnt += 1;
    xsm += x;
    x2sm += x*x;
  }

  public double getAve() {
    return (cnt > 0) ? xsm/cnt : 0;
  }

  public double getVar() {
    if (cnt <= 1) return 0;
    double av = xsm/cnt;
    return x2sm/cnt - av*av;
  }

  public double getStd() {
    return Math.sqrt( getVar() );
  }
}



