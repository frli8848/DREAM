//
// OpenCL delay-and-sum beamforming kernel
//

#define F_SFX(v) v##f // Suffix of numerical floating point float constants (= f suffix).

#define DAS_DATATYPE double
//#define DAS_DATATYPE float

__kernel void das_tfm(__global const double *Y, // Size: a_scan_len x num_t_elements*num_r_elements (=FMC)
                      const int a_scan_len,
                      __global const double *Gt, const int num_t_elements, // Size: num_t_elements x 3
                      __global const double *Gr, const int num_r_elements, // Size: num_r_elements x 3
                      __global const double *Ro, const int No,  // Size: No x 3
                      const double dt,
                      const double delay,
                      const double cp,
                      __global double *Im)
{
  int no = get_global_id(0); // Each thread computes one image point.
  const DAS_DATATYPE xo = (DAS_DATATYPE) Ro[no];
  const DAS_DATATYPE yo = (DAS_DATATYPE) Ro[no + No*1];
  const DAS_DATATYPE zo = (DAS_DATATYPE) Ro[no + No*2];

  const DAS_DATATYPE Fs_khz = (DAS_DATATYPE) (1.0/dt)*1000.0;
  const DAS_DATATYPE one_over_cp = (DAS_DATATYPE) 1.0/cp;
  const DAS_DATATYPE delay_ms = (DAS_DATATYPE) delay/1000.0;

  // Work on local data.
  //DAS_DATATYPE im = F_SFX(0.0); // For float
  DAS_DATATYPE im = 0.0; // For double

  __global double *y_p = (__global double *) Y;

  for (int n_t=0; n_t<num_t_elements; n_t++) {

    // Transmit
    DAS_DATATYPE gx_t = ((DAS_DATATYPE) Gt[n_t]) - xo;
    DAS_DATATYPE gy_t = ((DAS_DATATYPE) Gt[n_t + 1*num_t_elements]) - yo;
    DAS_DATATYPE gz_t = ((DAS_DATATYPE) Gt[n_t + 2*num_t_elements]) - zo;
    DAS_DATATYPE t_t = native_sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t) * one_over_cp; // [ms].
    t_t += delay_ms;            // Compensate for pulse system delay.

    for (int n_r=0; n_r<num_r_elements; n_r++) {

      // Recieve
      DAS_DATATYPE gx_r = ((DAS_DATATYPE) Gr[n_r]) - xo;
      DAS_DATATYPE gy_r = ((DAS_DATATYPE) Gr[n_r + 1*num_r_elements]) - yo;
      DAS_DATATYPE gz_r = ((DAS_DATATYPE) Gr[n_r + 2*num_r_elements]) - zo;
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
