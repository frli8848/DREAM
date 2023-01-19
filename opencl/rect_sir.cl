__kernel void rect_sir(__global const double *Ro,
                       int No,
                       double a_i,
                       double b_i,
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
  for (int it = 0; it < nt; it++) {
    h[it] = (double) 0.0;
  }

  int    it, k;
  double t, t_z;
  double pi;
  double tau_1, tau_2, tau_3, tau_4, a, b;
  double a_k=0, g_k=0, s_k=0, l_k=0;

  double cp_2pi = cp/(2.0*M_PI);

  // Convert to [m].
  xo = fabs(xo) / 1000.0;	// Can take abs due to symmetry.
  yo = fabs(yo) / 1000.0;	// Can take abs due to symmetry.
  zo = zo / 1000.0;
  a = a_i / 1000.0;
  b = b_i / 1000.0;

  for  (k=1; k<=4; k++) { // loop over all 4 sub-rectangles.

    if ( (xo <= a/2) && (yo <= b/2) ) {// Inside the aperture

      g_k = 1; // Add all.
      switch (k) {

      case 1:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 - yo);
        break;

      case 2:
        s_k = fabs(xo + a/2.0);
        l_k = fabs(b/2.0 - yo);
        break;

      case 3:
        s_k = fabs(xo + a/2.0);
        l_k = fabs(yo + b/2.0);
        break;

      case 4:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 - yo);
        break;

      default:
        break;
      }
    } // if in shadow.

    if ( (xo <= a/2.0) && (yo > b/2.0) ) {// Inside a/2 but outside b/2.

      switch (k) {

      case 1:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1.0; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      case 4:
        s_k = fabs(a/2.0 - xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1; // Inside => add.
        break;

      default:
        break;
      }
    } // if inside a/2 but outside b/2.

    if ( (xo > a/2.0) && (yo <= b/2.0) ) { // Inside b/2 but outside a/2.

      switch (k) {

      case 1:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(b/2.0 - yo);
        g_k = -1.0; // Outside => subtract.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo );
        l_k = fabs(b/2.0 - yo);
        g_k = 1.0; // Inside => add.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1.0; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(b/2.0 + yo);
        g_k = -1.0; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if inside b/2 but outside a/2.

    if ( (xo > a/2.0) && (yo > b/2.0) ) {// Outside both  a/2 and b/2.
      switch (k) {

      case 1:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(yo - b/2.0);
        g_k = 1.0; // Really outside but need to add since 1 is included both in 2 and 4.
        break;

      case 2:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(yo - b/2.0);
        g_k = -1.0; // Outside => subtract.
        break;

      case 3:
        s_k = fabs(a/2.0 + xo);
        l_k = fabs(b/2.0 + yo);
        g_k = 1.0; // Inside => add.
        break;

      case 4:
        s_k = fabs(xo - a/2.0);
        l_k = fabs(yo + b/2.0);
        g_k = -1.0; // Outside => subtract.
        break;

      default:
        break;
      }
    } // if outside both  a/2.0 and b/2.0.

    t_z = zo / cp;

    tau_1 = t_z;
    tau_2 = native_sqrt(zo*zo + s_k*s_k) / cp;
    tau_3 = native_sqrt(zo*zo + l_k*l_k) / cp;
    tau_4 = native_sqrt(zo*zo + l_k*l_k + s_k*s_k) / cp;

    for (it=0; it<nt; it++) {

      t = (((double) it) * dt + delay)/1.0e6; // in [s].

      a_k = 0.0;
      if ( (t >= tau_1) && (t <= tau_4) ) {
        a_k += cp/4.0;
      }

      if ( (t >= tau_2) && (t <= tau_4) ) {
        a_k -= cp_2pi * acos( s_k / (cp * native_sqrt(t*t - t_z*t_z)));
      }

      if ( (t >= tau_3) && (t <= tau_4) ) {
        a_k -= cp_2pi * acos( l_k / (cp * native_sqrt(t*t - t_z*t_z)));
      }

      if (a_k > 0.0) {
        h[it] += g_k * a_k; // This hits global memory
      }

    } // for it
  } // for k

  return;
}
