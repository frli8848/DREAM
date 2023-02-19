//
// OpenCL delay-and-sum synthetic aperture imaging beamforming kernels
// using uniform array element and observation point locations
// (on a rectangular/linear grids).
//

#define DAS_DATATYPE double

__kernel void das_uni_saft(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_elements
                           const int a_scan_len,
                           DAS_DATATYPE min_tr, DAS_DATATYPE pitch, DAS_DATATYPE max_tr,
                           DAS_DATATYPE min_Rx, DAS_DATATYPE dx, DAS_DATATYPE max_Rx,
                           DAS_DATATYPE min_Ry, DAS_DATATYPE dy, DAS_DATATYPE max_Ry,
                           DAS_DATATYPE min_Rz, DAS_DATATYPE dz, DAS_DATATYPE max_Rz,
                           const DAS_DATATYPE dt,
                           const DAS_DATATYPE delay,
                           const DAS_DATATYPE cp,
                           __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.

  // Array dims
  int num_elements = (int) ((max_tr - min_tr)/pitch+1.0);

  // Image dims
  int Nx = (int) ((max_Rx - min_Rx)/dx+1.0);
  int Ny = (int) ((max_Ry - min_Ry)/dy+1.0);
  int Nz = (int) ((max_Rz - min_Rz)/dz+1.0);

  // Convert linear indices to subscripts
  // (cf. the in2sub function in Octave/MATLAB).
  int y_quot = no / Ny;
  int ny = no - y_quot*Ny;

  int z_quot = y_quot / Nz;
  int nz = y_quot - z_quot*Nz;

  int x_quot = z_quot / Nx;
  int nx = z_quot - x_quot*Nx;

  DAS_DATATYPE xo = min_Rx + ((DAS_DATATYPE) nx)*dx;
  DAS_DATATYPE yo = min_Ry + ((DAS_DATATYPE) ny)*dy;
  DAS_DATATYPE zo = min_Rz + ((DAS_DATATYPE) nz)*dz;

  // Pre-compute this to avoid divisions in the inner loops.
  const DAS_DATATYPE Fs_khz = (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = 1.0/cp;
  const DAS_DATATYPE delay_ms = delay/1000.0;

  // Work on local data.
  DAS_DATATYPE im = 0.0;

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_tr=0; n_tr<num_elements; n_tr++) {

    // Transmit/Receive
    DAS_DATATYPE x_tr = min_tr + ((DAS_DATATYPE) n_tr)*pitch;
    DAS_DATATYPE gx_tr = x_tr - xo;
    //DAS_DATATYPE gy_tr = yo; // Assume y_tr = 0.0;
    //DAS_DATATYPE gz_tr = zo; // Assume z_tr = 0.0;
    DAS_DATATYPE t_tr = native_sqrt(gx_tr*gx_tr + yo*yo + zo*zo) * one_over_cp; // [ms].
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

__kernel void das_uni_tfm(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                          const int a_scan_len,
                          DAS_DATATYPE min_t, DAS_DATATYPE pitch_t, DAS_DATATYPE max_t,
                          DAS_DATATYPE min_r, DAS_DATATYPE pitch_r, DAS_DATATYPE max_r,
                          DAS_DATATYPE min_Rx, DAS_DATATYPE dx, DAS_DATATYPE max_Rx,
                          DAS_DATATYPE min_Ry, DAS_DATATYPE dy, DAS_DATATYPE max_Ry,
                          DAS_DATATYPE min_Rz, DAS_DATATYPE dz, DAS_DATATYPE max_Rz,
                          const DAS_DATATYPE dt,
                          const DAS_DATATYPE delay,
                          const DAS_DATATYPE cp,
                          __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.

  //
  // Get coordinates for the obervarion (image) point, ro, at
  // at linear index no.
  //

  // Array dims
  int num_t_elements = (int) ((max_t - min_t)/pitch_t+1.0);
  int num_r_elements = (int) ((max_r - min_r)/pitch_r+1.0);

  // Image dims
  int Nx = (int) ((max_Rx - min_Rx)/dx+1.0);
  int Ny = (int) ((max_Ry - min_Ry)/dy+1.0);
  int Nz = (int) ((max_Rz - min_Rz)/dz+1.0);

  // Convert linear indices to subscripts
  // (cf. the in2sub function in Octave/MATLAB).
  int y_quot = no / Ny;
  int ny = no - y_quot*Ny;

  int z_quot = y_quot / Nz;
  int nz = y_quot - z_quot*Nz;

  int x_quot = z_quot / Nx;
  int nx = z_quot - x_quot*Nx;

  DAS_DATATYPE xo = min_Rx + ((DAS_DATATYPE) nx)*dx;
  DAS_DATATYPE yo = min_Ry + ((DAS_DATATYPE) ny)*dy;
  DAS_DATATYPE zo = min_Rz + ((DAS_DATATYPE) nz)*dz;

  // Pre-compute this to avoid divisions in the inner loops.
  const DAS_DATATYPE Fs_khz = (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = 1.0/cp;
  const DAS_DATATYPE delay_ms = delay/1000.0;

  // Work on local data.
  DAS_DATATYPE im = 0.0;

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    DAS_DATATYPE xt = min_t + ((DAS_DATATYPE) n_t)*pitch_t;
    DAS_DATATYPE gx_t = xt - xo;
    //DAS_DATATYPE gy_t = yo; // Assume yt = 0.0;
    //DAS_DATATYPE gz_t = zo; // Assume zt = 0.0;
    //DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp; // [ms].
    DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + yo*yo + zo*zo) * one_over_cp; // [ms].
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      DAS_DATATYPE xr = min_r + ((DAS_DATATYPE) n_r)*pitch_r;
      DAS_DATATYPE gx_r = xr - xo;
      //DAS_DATATYPE gy_r = yo;  // Assume yr = 0.0;
      //DAS_DATATYPE gz_r = zo;  // Assume zt = 0.0;
      //DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r) * one_over_cp;
      DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + yo*yo + zo*zo) * one_over_cp; // [ms].

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

__kernel void das_uni_rca(__global const DAS_DATATYPE *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                          const int a_scan_len,
                          DAS_DATATYPE min_t, DAS_DATATYPE pitch_t, DAS_DATATYPE max_t,
                          DAS_DATATYPE min_r, DAS_DATATYPE pitch_r, DAS_DATATYPE max_r,
                          DAS_DATATYPE min_Rx, DAS_DATATYPE dx, DAS_DATATYPE max_Rx,
                          DAS_DATATYPE min_Ry, DAS_DATATYPE dy, DAS_DATATYPE max_Ry,
                          DAS_DATATYPE min_Rz, DAS_DATATYPE dz, DAS_DATATYPE max_Rz,
                          const DAS_DATATYPE dt,
                          const DAS_DATATYPE delay,
                          const DAS_DATATYPE cp,
                          __global DAS_DATATYPE *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.

  //
  // Get coordinates for the obervarion (image) point, ro, at
  // at linear index no.
  //

  // Array dims
  int num_t_elements = (int) ((max_t - min_t)/pitch_t+1.0);
  int num_r_elements = (int) ((max_r - min_r)/pitch_r+1.0);

  // Image dims
  int Nx = (int) ((max_Rx - min_Rx)/dx+1.0);
  int Ny = (int) ((max_Ry - min_Ry)/dy+1.0);
  int Nz = (int) ((max_Rz - min_Rz)/dz+1.0);

  // Convert linear indices to subscripts
  // (cf. the in2sub function in Octave/MATLAB).
  int y_quot = no / Ny;
  int ny = no - y_quot*Ny;

  int z_quot = y_quot / Nz;
  int nz = y_quot - z_quot*Nz;

  int x_quot = z_quot / Nx;
  int nx = z_quot - x_quot*Nx;

  DAS_DATATYPE xo = min_Rx + ((DAS_DATATYPE) nx)*dx;
  DAS_DATATYPE yo = min_Ry + ((DAS_DATATYPE) ny)*dy;
  DAS_DATATYPE zo = min_Rz + ((DAS_DATATYPE) nz)*dz;

  // Element (stripe) lengths
  const DAS_DATATYPE gt_y_min = min_r;
  const DAS_DATATYPE gt_y_max = max_r;
  const DAS_DATATYPE gr_x_min = min_t;
  const DAS_DATATYPE gr_x_max = max_t;

  // Pre-compute this to avoid divisions in the inner loops.
  const DAS_DATATYPE Fs_khz = (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = 1.0/cp;
  const DAS_DATATYPE delay_ms = delay/1000.0;

  // Work on local data.
  DAS_DATATYPE im = 0.0;

  __global DAS_DATATYPE *y_p = (__global DAS_DATATYPE *) Y;

  for (int n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    DAS_DATATYPE xt = min_t + ((DAS_DATATYPE) n_t)*pitch_t;
    DAS_DATATYPE gx_t = xt - xo;
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

    //DAS_DATATYPE gz_t = zo; // Assume zt = 0.0;
    DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + gy_t*gy_t + zo*zo) * one_over_cp; // [ms].
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      DAS_DATATYPE gx_r = 0.0;
      if ( (xo < gr_x_min) || xo > gr_x_max) {
        // We are outside the stripe aperture and
        // we use the distance to the edge of the stripe.
        if (xo < gr_x_min) {
          gx_r = gr_x_min - xo;
        } else { // yo > gt_y_max
          gx_r = gr_x_max - xo;
        }
      }
      DAS_DATATYPE yr = min_r + ((DAS_DATATYPE) n_r)*pitch_r;
      DAS_DATATYPE gy_r = yr - yo;
      //DAS_DATATYPE gz_r = zo; // Assume zr = 0.0;
      DAS_DATATYPE t_r = native_sqrt(gx_r*gx_r + gy_r*gy_r + zo*zo) * one_over_cp;

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
