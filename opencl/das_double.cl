//
// OpenCL delay-and-sum synthetic aperture imaging beamforming kernels
// using arbitrary array element and observation point locations.
//

#define DAS_DATATYPE double

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

  const DAS_DATATYPE Fs_khz = (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = 1.0/cp;
  const DAS_DATATYPE delay_ms = delay/1000.0;

  // Work on local data.
  DAS_DATATYPE im = 0.0;

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

  const DAS_DATATYPE Fs_khz = (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = 1.0/cp;
  const DAS_DATATYPE delay_ms = delay/1000.0;

  // Work on local data.
  DAS_DATATYPE im = 0.0; // For double

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

__kernel void das_rca_coltx(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                            const int a_scan_len,
                            __global const DAS_DATATYPE *G_col, const int num_cols, // Size: num_cols x 3
                            __global const DAS_DATATYPE *G_row, const int num_rows, // Size: num_rows x 3
                            __global const DAS_DATATYPE *Ro, const int No,  // Size: No x 3
                            const DAS_DATATYPE dt,
                            const DAS_DATATYPE delay,
                            const DAS_DATATYPE cp,
                            __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = (DAS_DATATYPE) Ro[no];
  const DAS_DATATYPE yo = (DAS_DATATYPE) Ro[no + No*1];
  const DAS_DATATYPE zo = (DAS_DATATYPE) Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (DAS_DATATYPE) (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = (DAS_DATATYPE) 1.0/cp;
  const DAS_DATATYPE delay_ms = (DAS_DATATYPE) delay/1000.0;

  // Element (stripe) lengths
  const DAS_DATATYPE gc_y_min = G_row[0 + 1*num_cols];
  const DAS_DATATYPE gc_y_max = G_row[num_cols-1 + 1*num_cols];
  const DAS_DATATYPE gr_x_min = G_col[0];
  const DAS_DATATYPE gr_x_max = G_col[num_rows-1];

  // Work on local data.
  DAS_DATATYPE im = 0.0;

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_c=0; n_c<num_cols; n_c++) {

    // Cols
    DAS_DATATYPE gx_c = ((DAS_DATATYPE) G_col[n_c]) - xo;
    DAS_DATATYPE gy_c;
    if ( (yo >= gc_y_min) && yo <= gc_y_max) {
      gy_c = 0.0;               // We are inside the stripe aperture.
    } else {    // Here we use the distance to the edge of the stripe.
      if (yo < gc_y_min) {
        gy_c = gc_y_min - yo;
      } else { // yo > gc_y_max
        gy_c = gc_y_max - yo;
      }
    }
    DAS_DATATYPE gz_c = G_col[n_c + 2*num_cols] - zo;
    DAS_DATATYPE t_c = native_sqrt(gx_c*gx_c + gy_c*gy_c + gz_c*gz_c) * one_over_cp; // [ms].
    t_c += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_rows; n_r++) {

      // Rows
      DAS_DATATYPE gx_r = 0.0;
      if ( (xo < gr_x_min) || xo > gr_x_max) {
        // We are outside the stripe aperture and
        // we use the distance to the edge of the stripe.
        if (xo < gr_x_min) {
          gx_r = gr_x_min - xo;
        } else { // yo > gc_y_max
          gx_r = gr_x_max - xo;
        }
      }
      DAS_DATATYPE gy_r = G_row[n_r + 1*num_rows] - yo;
      DAS_DATATYPE gz_r = G_row[n_r + 2*num_rows] - zo;
      DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

      DAS_DATATYPE t_dp = t_c + t_r; // Double-path travel time.
      int k = (int) (t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += y_p[k];
      }

      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  // Just write to global memory once!
  Im[no] = im;
}

__kernel void das_rca_rowtx(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_cols*num_rows (=FMC)
                            const int a_scan_len,
                            __global const DAS_DATATYPE *G_col, const int num_cols, // Size: num_cols x 3
                            __global const DAS_DATATYPE *G_row, const int num_rows, // Size: num_rows x 3
                            __global const DAS_DATATYPE *Ro, const int No,  // Size: No x 3
                            const DAS_DATATYPE dt,
                            const DAS_DATATYPE delay,
                            const DAS_DATATYPE cp,
                            __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = (DAS_DATATYPE) Ro[no];
  const DAS_DATATYPE yo = (DAS_DATATYPE) Ro[no + No*1];
  const DAS_DATATYPE zo = (DAS_DATATYPE) Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (DAS_DATATYPE) (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = (DAS_DATATYPE) 1.0/cp;
  const DAS_DATATYPE delay_ms = (DAS_DATATYPE) delay/1000.0;

  // Element (stripe) lengths
  const DAS_DATATYPE gc_y_min = G_row[0 + 1*num_cols];
  const DAS_DATATYPE gc_y_max = G_row[num_cols-1 + 1*num_cols];
  const DAS_DATATYPE gr_x_min = G_col[0];
  const DAS_DATATYPE gr_x_max = G_col[num_rows-1];

  // Work on local data.
  DAS_DATATYPE im = 0.0;

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_r=0; n_r<num_rows; n_r++) {

    // Rows
    DAS_DATATYPE gx_r = 0.0;
    if ( (xo < gr_x_min) || xo > gr_x_max) {
      // We are outside the stripe aperture and
      // we use the distance to the edge of the stripe.
      if (xo < gr_x_min) {
        gx_r = gr_x_min - xo;
      } else { // yo > gc_y_max
        gx_r = gr_x_max - xo;
      }
    }
    DAS_DATATYPE gy_r = G_row[n_r + 1*num_rows] - yo;
    DAS_DATATYPE gz_r = G_row[n_r + 2*num_rows] - zo;
    DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;

    for (int n_c=0; n_c<num_cols; n_c++) {

      // Cols
      DAS_DATATYPE gx_c = ((DAS_DATATYPE) G_col[n_c]) - xo;
      DAS_DATATYPE gy_c;
      if ( (yo >= gc_y_min) && yo <= gc_y_max) {
        gy_c = 0.0;               // We are inside the stripe aperture.
      } else {    // Here we use the distance to the edge of the stripe.
        if (yo < gc_y_min) {
          gy_c = gc_y_min - yo;
        } else { // yo > gc_y_max
          gy_c = gc_y_max - yo;
        }
      }
      DAS_DATATYPE gz_c = G_col[n_c + 2*num_cols] - zo;
      DAS_DATATYPE t_c = native_sqrt(gx_c*gx_c + gy_c*gy_c + gz_c*gz_c) * one_over_cp; // [ms].
      t_c += delay_ms;            // Compensate for pulse system delay.

      DAS_DATATYPE t_dp = t_c + t_r; // Double-path travel time.
      int k = (int) (t_dp*Fs_khz);

      if ((k < a_scan_len) && (k >= 0)) {
        im += y_p[k];
      }

      y_p += a_scan_len; // Jump to the next A-scan
    }
  }

  // Just write to global memory once!
  Im[no] = im;
}
