__kernel void dreamrect(__global const double *Ro,
                        int No,
                        double a,
                        double b,
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

  double xsmin = -a/2.0;
  double xsmax =  a/2.0;
  double ysmin = -b/2.0;
  double ysmax =  b/2.0;

  double ds = dx * dy;

  double rz = zo;
  double ys = ysmin + dy / 2.0;

  while (ys <= ysmax) {

    double ry = yo - ys;
    double xs = xsmin + dx / 2.0;

    while (xs <= xsmax) {

      double rx = xo - xs;
      double r = sqrt(rx*rx + ry*ry + rz*rz);

      double ai = v * ds / (2.0*M_PI * r);
      ai /= dt;
      ai *= 1000.0;		// Convert to SI units.

      double t = r * 1000.0/cp;	// Propagation delay in micro seconds.
      int it = (int) rint((t - delay)/dt); // Sample index.

      // Check if index is out of bounds.
      if ( (it < nt) && (it >= 0) ) {
        h[it] += ai; // TODO here we hit global memory - try to avoid this in the inner loop.
      }

      xs += dx;
    }
    ys += dy;
  }
}
