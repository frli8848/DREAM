__kernel void dreamcirc(__global const double *Ro,
                        int No,
                        double r,
                        double dx,
                        double dy,
                        double dt,
                        int nt,
                        double delay,
                        double v,
                        double cp,
                        __global double *H)
{
  int i, it;
  double t;
  double ai;
  double ri, x, y;
  double rx,ry,rz;

  int no = get_global_id(0);
  double xo = Ro[no];
  double yo = Ro[no + No*1];
  double zo = Ro[no + No*2];

  double pi = 4.0 * atan(1.0);
  double ds = dx * dy;

  rs = sqrt(r*r - (ys-y)*(ys-y));
  xsmin = -rs + xs;
  xsmax = rs + xs;

  ry = yo - y;

  x = xsmin + dx / 2.0;

  __global double *h = &H[0+nt*no];
  for (i = 0; i < nt; i++) {
    h[i] = 0.0;
  }

  y = ysmin + dy/2;
  while (y <= ysmax) {

    //xlimit_circ(y, r, xs, ys, &xsmin, &xsmax);
    rs = sqrt(r*r - (ys-y)*(ys-y));
    xsmin = -rs + xs;
    xsmax = rs + xs;

    ry = yo - y;

    x = xsmin + dx / 2.0;
    while (x <= xsmax) {

      //distance(xo, yo, zo, x, y, &ri);
      rx = xo - x;
      //ry = yo - y; // Moved outside this loop.
      //rz = zo;
      ri = sqrt(rx*rx + ry*ry + zo*zo);

      ai = v * ds / (2*pi * ri);
      ai /= dt;
      ai *= 1000; // Convert to SI units.

      // Propagation delay in micro seconds.
      t = ri * 1000.0/cp;
      it = (int) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai;
      }

      x += dx;
    }

    y += dy;
  }
}
