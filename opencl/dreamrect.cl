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
  int i, it;
  double t;
  double ai;
  double ri, x, y;
  double rx,ry,rz;

  double xsmin = -a/2.0;
  double xsmax =  a/2.0;
  double ysmin = -b/2.0;
  double ysmax =  b/2.0;

  double pi = atan( (double) 1.0) * 4.0;
  double ds = dx * dy;

  int no = get_global_id(0);
  double xo = Ro[no];
  double yo = Ro[no + No*1];
  double zo = Ro[no + No*2];

  __global double *h = &H[0+nt*no];
  for (i = 0; i < nt; i++) {
    h[i] = 0.0;
  }

  rz = zo;
  y = ysmin + dy / 2.0;

  while (y <= ysmax) {
    ry = yo - y;
    x = xsmin + dx / 2.0;

    while (x <= xsmax) {

      rx = xo - x;
      ri = sqrt(rx*rx + ry*ry + rz*rz);

      ai = v * ds / (2*pi * ri);
      ai /= dt;
      ai *= 1000.0;		// Convert to SI units.

      t = ri * 1000.0/cp;	// Propagation delay in micro seconds.
      it = (int) rint((t - delay)/dt); // Sample index.

      // Check if index is out of bounds.
      if ( (it < nt) && (it >= 0) ) {
        h[it] += ai; // TODO here we hit global memory - try to avoid this in the inner loop.
      }

      x += dx;
    }
    y += dy;
  }

}
