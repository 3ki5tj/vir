class Matrix3D {
  double xx, xy, xz, xo;
  double yx, yy, yz, yo;
  double zx, zy, zz, zo;
  static final double pi = 3.14159265;

  /** Create a new unit matrix */
  Matrix3D() {
    xx = 1.0;
    yy = 1.0;
    zz = 1.0;
  }

  /** Scale along each axis independently */
  void scale(double xf, double yf, double zf) {
    xx *= xf; xy *= xf; xz *= xf; xo *= xf;
    yx *= yf; yy *= yf; yz *= yf; yo *= yf;
    zx *= zf; zy *= zf; zz *= zf; zo *= zf;
  }

  /** Translate the origin */
  void translate(double x, double y, double z) {
    xo += x;
    yo += y;
    zo += z;
  }

  /** Rotate theta degrees around the y axis */
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

  /** Rotate theta degrees about the x axis */
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

  /** Rotate theta degrees about the z axis */
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

  /** Recover the unit matrix */
  void unit() {
    xo = 0; xx = 1; xy = 0; xz = 0;
    yo = 0; yx = 0; yy = 1; yz = 0;
    zo = 0; zx = 0; zy = 0; zz = 1;
  }

  /** Transform np points from v into tv.
   *  v contains the input coordinates in floating point.
   *  Three successive entries in the array constitute a point.
   *  tv ends up holding the transformed points as integers;
   *  three successive entries per point */
  void transform(double v[][], int tv[][], int np) {
    // np can be different from v.length
    for (int i = 0; i < np; i++) {
      double x = v[i][0], y = v[i][1], z = v[i][2];
      tv[i][0] = (int) (xx * x + xy * y + xz * z + xo);
      tv[i][1] = (int) (yx * x + yy * y + yz * z + yo);
      tv[i][2] = (int) (zx * x + zy * y + zz * z + zo);
    }
  }
}




