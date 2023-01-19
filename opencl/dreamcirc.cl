__kernel void dreamcirc(__global const double *Ro,
                        int No,
                        double R,
                        double dx,
                        double dy,
                        double dt,
                        int nt,
                        double delay,
                        double v,
                        double cp,
                        __global double *H)
{
  int no = get_global_id(0);
  double xo = Ro[no];
  double yo = Ro[no + No*1];
  double zo = Ro[no + No*2];

  __global double *h = &H[0+nt*no];
  for (int i = 0; i < nt; i++) {
    h[i] = 0.0;
  }

  double ds = dx * dy;

  // y-dim integration bounds
  double ysmin = -R;
  double ysmax =  R;

  double ys = ysmin + dy/2.0;
  while (ys <= ysmax) {

    double rxs = sqrt(R*R - ys*ys);
    double xsmin = -rxs;
    double xsmax = rxs;

    double ry = yo - ys;

    double xs = xsmin + dx/2.0;
    while (xs <= xsmax) {

      // Compute the distance (length) from an observation point (xo,yo,zo)
      // to a point (xs,ys) on the transducer surface.
      double rx = xo - xs;
      double r = sqrt(rx*rx + ry*ry + zo*zo);

      double ai = v * ds / (2*M_PI * r);
      ai /= dt;
      ai *= 1000.0; // Convert to SI units.

      // Propagation delay in micro seconds.
      double t = r * 1000.0/cp;
      int it = (int) rint((t - delay)/dt);

      // Check if index is out of bounds.
      if ((it < nt) && (it >= 0)) {
        h[it] += ai; // TODO here we hit global memory - try to avoid this in the inner loop.
      }

      xs += dx;
    }
    ys += dy;
  }
}
