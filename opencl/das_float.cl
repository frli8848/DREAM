//
// OpenCL delay-and-sum synthetic aperture imaging beamforming kernels
// using arbitrary array element and observation point locations.
//

#define F_SFX(v) v##f // Suffix of numerical floating point float constants (= f suffix).
#define DAS_DATATYPE float

__kernel void das_saft(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_elements
                       const int a_scan_len,
                       __global const DAS_DATATYPE *G, const int num_elements, // Size: num_elements x 3
                       __global const DAS_DATATYPE *Ro, const int No,  // Size: No x 3
                       const DAS_DATATYPE dt,
                       const DAS_DATATYPE delay,
                       const DAS_DATATYPE cp,
                       __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = Ro[no];
  const DAS_DATATYPE yo = Ro[no + No*1];
  const DAS_DATATYPE zo = Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (F_SFX(1.0)/dt)*F_SFX(1000.0);
  const DAS_DATATYPE one_over_cp = F_SFX(1.0)/cp;
  const DAS_DATATYPE delay_ms = delay/F_SFX(1000.0);

  // Work on local data.
  DAS_DATATYPE im = F_SFX(0.0); // For float

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_tr=0; n_tr<num_elements; n_tr++) {

    // Transmit/Receive
    DAS_DATATYPE gx_tr = G[n_tr] - xo;
    DAS_DATATYPE gy_tr = G[n_tr + 1*num_elements] - yo;
    DAS_DATATYPE gz_tr = G[n_tr + 2*num_elements] - zo;
    DAS_DATATYPE t_tr = native_sqrt(gx_tr*gx_tr + gy_tr*gy_tr + gz_tr*gz_tr) * one_over_cp; // [ms].
    t_tr += delay_ms;           // Compensate for pulse system delay.

    DAS_DATATYPE t_dp = 2.0*t_tr; // Double-path travel time.
    int k = (int) (t_dp*Fs_khz);

    if ((k < a_scan_len) && (k >= 0)) {
      im += (DAS_DATATYPE) y_p[k];
    }

    y_p += a_scan_len; // Jump to the next A-scan
  }

  // Just write to global memory once!
  Im[no] = im;
}

__kernel void das_tfm(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                      const int a_scan_len,
                      __global const DAS_DATATYPE *Gt, const int num_t_elements, // Size: num_t_elements x 3
                      __global const DAS_DATATYPE *Gr, const int num_r_elements, // Size: num_r_elements x 3
                      __global const DAS_DATATYPE *Ro, const int No,  // Size: No x 3
                      const DAS_DATATYPE dt,
                      const DAS_DATATYPE delay,
                      const DAS_DATATYPE cp,
                      __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = Ro[no];
  const DAS_DATATYPE yo = Ro[no + No*1];
  const DAS_DATATYPE zo = Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (F_SFX(1.0)/dt)*F_SFX(1000.0);
  const DAS_DATATYPE one_over_cp = F_SFX(1.0)/cp;
  const DAS_DATATYPE delay_ms = delay/F_SFX(1000.0);

  // Work on local data.
  DAS_DATATYPE im = F_SFX(0.0); // For float

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    DAS_DATATYPE gx_t = Gt[n_t] - xo;
    DAS_DATATYPE gy_t = Gt[n_t + 1*num_t_elements] - yo;
    DAS_DATATYPE gz_t = Gt[n_t + 2*num_t_elements] - zo;
    DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp; // [ms].
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      DAS_DATATYPE gx_r = Gr[n_r] - xo;
      DAS_DATATYPE gy_r = Gr[n_r + 1*num_r_elements] - yo;
      DAS_DATATYPE gz_r = Gr[n_r + 2*num_r_elements] - zo;
      DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      DAS_DATATYPE t_dp = t_t + t_r; // Double-path travel time.
      int k = (int) (t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += (DAS_DATATYPE) y_p[k];
      }

      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  // Just write to global memory once!
  Im[no] = im;
}

__kernel void das_rca(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                      const int a_scan_len,
                      __global const DAS_DATATYPE *Gt, const int num_t_elements, // Size: num_t_elements x 3
                      __global const DAS_DATATYPE *Gr, const int num_r_elements, // Size: num_r_elements x 3
                      __global const DAS_DATATYPE *Ro, const int No,  // Size: No x 3
                      const DAS_DATATYPE dt,
                      const DAS_DATATYPE delay,
                      const DAS_DATATYPE cp,
                      __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = Ro[no];
  const DAS_DATATYPE yo = Ro[no + No*1];
  const DAS_DATATYPE zo = Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (DAS_DATATYPE) (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = (DAS_DATATYPE) 1.0/cp;
  const DAS_DATATYPE delay_ms = (DAS_DATATYPE) delay/1000.0;

  // Element (stripe) lengths
  const DAS_DATATYPE gt_y_min = Gr[0 + 1*num_t_elements];
  const DAS_DATATYPE gt_y_max = Gr[num_t_elements-1 + 1*num_t_elements];
  const DAS_DATATYPE gr_x_min = Gt[0];
  const DAS_DATATYPE gr_x_max = Gt[num_r_elements-1];

  // Work on local data.
  DAS_DATATYPE im = F_SFX(0.0); // For float

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    DAS_DATATYPE gx_t = Gt[n_t] - xo;
    DAS_DATATYPE gy_t;
    if ( (yo >= gt_y_min) && yo <= gt_y_max) {
      gy_t = 0.0;               // We are inside the stripe aperture.
    } else {    // Here we use the distance to the edge of the stripe.
      if (yo < gt_y_min) {
        gy_t = gt_y_min - yo;
      } else { // yo > gt_y_max
        gy_t = gt_y_max - yo;
      }
    }
    DAS_DATATYPE gz_t = Gt[n_t + 2*num_t_elements] - zo;
    DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp; // [ms].
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      DAS_DATATYPE gx_r = F_SFX(0.0);
      if ( (xo < gr_x_min) || xo > gr_x_max) {
        // We are outside the stripe aperture and
        // we use the distance to the edge of the stripe.
        if (xo < gr_x_min) {
          gx_r = gr_x_min - xo;
        } else { // yo > gt_y_max
          gx_r = gr_x_max - xo;
        }
      }
      DAS_DATATYPE gy_r = Gr[n_r + 1*num_r_elements] - yo;
      DAS_DATATYPE gz_r = Gr[n_r + 2*num_r_elements] - zo;
      DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      DAS_DATATYPE t_dp = t_t + t_r; // Double-path travel time.
      int k = (int) (t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += (DAS_DATATYPE) y_p[k];
      }

      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  // Just write to global memory once!
  Im[no] = im;
}
