/***
*
* Copyright (C) 2023 Fredrik Lingvall
*
* This file is part of the DREAM Toolbox.
*
* The DREAM Toolbox is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by the
* Free Software Foundation; either version 2, or (at your option) any
* later version.
*
* The DREAM Toolbox is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
* for more details.
*
* You should have received a copy of the GNU General Public License
* along with the DREAM Toolbox; see the file COPYING.  If not, write to the
* Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
* 02110-1301, USA.
*
***/

#include<csignal>

#include "das_uni.h"
#include "arg_parser.h"

#include <octave/oct.h>

// Persistent smart pointer to DAS_UNI object.
// This one is only cleared when we do a
// octave:1> clear das_uni
std::unique_ptr<DAS_UNI<double>> das_uni_d=nullptr;
std::unique_ptr<DAS_UNI<float>> das_uni_f=nullptr;

/***
 *
 * das_uni - Octave (oct) gateway function for das_uni (delay-and-sum using uniform grids).
 *
 ***/

DEFUN_DLD (das_uni, args, nlhs,
           "-*- texinfo -*-\n\
@deftypefn {Loadable Function} {}  [Im] = das_uni(Y,gt,gr,Ro,dt,delay,cp,method).\n \
\n\
DAS_UNI Computes the delay-and-sum processed reconstruction (beamformed image) for\n\
three different array geometries: SAFT, TFM, and RCA using uniform grids both for the\n\
array element positions and for the image points. In the synthtic aperture\n\
focusing techinque (SAFT) one use one transmitter and one reciever (the same)\n\
and moves that along the array aperture. In the total focusing technique (TFM) one\n\
first transmit with the first element and then recieve with all elements, transmit with the \n\
second element and again receive with all elements and so on until the last transmit element (ie.,\n\
transmit with all - receive with all). The last method, row-column adressed (RCA) array, is \n\
a variant of TFM where the elements are arranged in a crossed layout to form a 2D array.\n\
\n\
Data matrix:\n\
\n\
@table @code\n\
@item Y\n\
A K x N matrix where K is the A-scan length and the number of A-scans, N, depends\n\
on DAS algorithm selected (see below).\n\
@end table\n\
\n\
Transmit element vector [mm]:\n\
\n\
@table @code\n\
\n\
@item gt\n\
A 3 element vector, gt = [min_t, pitch_t, max_t].\n\
@end table\n\
\n\
Receive element vector [mm]:\n\
\n\
@table @code\n\
@item gr\n\
A 3 element vector, gr = [min_r, pitch_r, max_r].\n\
@end table\n\
\n\
Observation (image) point matrix [mm]:\n\
\n\
@table @code\n\
@item Ro\n\
A 3 x 3 matrix, Ro = [min_Rx, dx, max_Rx; min_Ry, dy, max_Ry; min_Rz, dz, max_Rz].\n\
@end table\n\
\n\
Sampling parameter dt: \n\
\n\
@table @code\n\
@item dt\n\
Temporal discretization period (= 1/sampling freq) [us].\n\
@end table\n\
\n\
Data start and pulse delay compensation:\n\
\n\
@table @code\n\
@item  delay\n\
Scalar delay for all observation points or a vector with individual delays for each observation point [us].\n\
@end table\n\
\n\
Sound speed:\n\
\n\
@table @code\n\
@item cp\n\
Sound velocity of the medium [m/s].\n\
\n\
@end table\n\
DAS algorithm:\n\
das_met is a text string parameter for selecting DAS algorithm, options are:\n\
\n\
@table @code\n\
@item 'saft'\n\
When SAFT is selected data Y must be an K x Lt (SAFT is also selected if Gr=[]).\n\
@item 'tfm'\n\
When TFM is selected a linear array is assumed and data Y must be a \n\
K x Lt*Lr matrix.\n\
@item 'rca_coltx' or 'rca_coltx'\n\
When RCA is selected a 2D RCA array is assumed and data Y must be a \n\
K x Lt*Lr matrix.\n\
@end table\n\
\n\
das_uni is an oct-function that is a part of the DREAM Toolbox available at\n\
@url{https://github.com/frli8848/DREAM}.\n\
\n\
Copyright @copyright{} 2023 Fredrik Lingvall.\n\
@end deftypefn")
{
  octave_value_list oct_retval;

  int nrhs = args.length ();

  ArgParser ap;

  //
  // Check for proper number of arguments
  //

  if (nrhs != 8) {
    error("das_uni requires 8 input arguments!");
    return oct_retval;
  }

  if (nlhs > 1) {
    error("Too many output arguments for das_uni!");
    return oct_retval;
  }

  //
  // Check if we are using single or double precision processsing
  //

  // Args 1 to 7 must have the same datatype (float or double).
  bool use_float = args(0).is_single_type();
  if (use_float) {
    if (args(1).is_single_type() != use_float ||
        args(2).is_single_type() != use_float ||
        args(3).is_single_type() != use_float ||
        args(4).is_single_type() != use_float ||
        args(5).is_single_type() != use_float ||
        args(6).is_single_type() != use_float) {
      error("First arg is single precision but one of args 2 to 7 is not!");
      return oct_retval;
    }
  } else {
    if (!(args(1).is_double_type() &&
          args(2).is_double_type() &&
          args(3).is_double_type() &&
          args(4).is_double_type() &&
          args(5).is_double_type() &&
          args(6).is_double_type()) ) {
      error("First arg is double precision but one of args 2 to 7 is not!");
      return oct_retval;
    }
  }

  //
  // Data
  //

  dream_idx_type a_scan_len = mxGetM(0); // A-scan length
  dream_idx_type num_a_scans = mxGetN(0);

  const float *Yf = nullptr;
  const double *Yd = nullptr;

  if (use_float) {
    const FloatMatrix tmp0f= args(0).float_matrix_value();
    Yf = (float*) tmp0f.data();
  } else {
    const Matrix tmp0d = args(0).matrix_value();
    Yd = (double*) tmp0d.data();
  }

  //
  // Transmit array
  //

  if (mxGetM(1)*mxGetN(1) != 3) {
    error("Argument 2 must be a 3 element vector");
    return oct_retval;
  }

  float min_t_f = 0.0;
  float pitch_t_f = 0.0;
  float max_t_f = 0.0;

  double min_t_d = 0.0;
  double pitch_t_d = 0.0;
  double max_t_d = 0.0;

  if (use_float) {
    const FloatMatrix tmp1f = args(1).float_matrix_value();
    min_t_f = tmp1f.data()[0];
    pitch_t_f = tmp1f.data()[1];
    max_t_f = tmp1f.data()[2];
  } else {
    const Matrix tmp1d = args(1).matrix_value();
    min_t_d = tmp1d.data()[0];
    pitch_t_d = tmp1d.data()[1];
    max_t_d = tmp1d.data()[2];
  }

  //
  // Recieve array
  //

  float min_r_f = 0.0;
  float pitch_r_f = 0.0;
  float max_r_f = 0.0;

  double min_r_d = 0.0;
  double pitch_r_d = 0.0;
  double max_r_d = 0.0;

  if (mxGetM(2)*mxGetN(2) != 0) { // SAFT is using an empty gr.

    if (mxGetM(2)*mxGetN(2) != 3) {
      error("Argument 3 must be a 3 element vector!");
      return oct_retval;
    }

    if (use_float) {
      const FloatMatrix tmp2f = args(2).float_matrix_value();
      min_r_f = tmp2f.data()[0];
      pitch_r_f = tmp2f.data()[1];
      max_r_f = tmp2f.data()[2];
    } else {
      const Matrix tmp2d = args(2).matrix_value();
      min_r_d = tmp2d.data()[0];
      pitch_r_d = tmp2d.data()[1];
      max_r_d = tmp2d.data()[2];
    }

  }

  //
  // Observation points (in a uniform grid).
  //

  // Check that arg is a 3 x 3 matrix.
  if (mxGetM(3) != 3 || mxGetN(3) != 3) {
    error("Argument 4 must be a 3 x 3 matrix!");
    return oct_retval;
  }

  dream_idx_type Nx = 0, Ny = 0, Nz = 0;

  float min_Rx_f=0.0, dx_f=0.0, max_Rx_f=0.0;
  float min_Ry_f=0.0, dy_f=0.0, max_Ry_f=0.0;
  float min_Rz_f=0.0, dz_f=0.0, max_Rz_f=0.0;

  double min_Rx_d=0.0, dx_d=0.0, max_Rx_d=0.0;
  double min_Ry_d=0.0, dy_d=0.0, max_Ry_d=0.0;
  double min_Rz_d=0.0, dz_d=0.0, max_Rz_d=0.0;

  if (use_float) {
    const FloatMatrix tmp3f= args(3).float_matrix_value();
    min_Rx_f=tmp3f.data()[0]; dx_f=tmp3f.data()[3]; max_Rx_f=tmp3f.data()[6];
    min_Ry_f=tmp3f.data()[1]; dy_f=tmp3f.data()[4]; max_Ry_f=tmp3f.data()[7];
    min_Rz_f=tmp3f.data()[2]; dz_f=tmp3f.data()[5]; max_Rz_f=tmp3f.data()[8];
    Nx = (dream_idx_type) ((max_Rx_f - min_Rx_f)/dx_f+1.0);
    Ny = (dream_idx_type) ((max_Ry_f - min_Ry_f)/dy_f+1.0);
    Nz = (dream_idx_type) ((max_Rz_f - min_Rz_f)/dz_f+1.0);
  } else {
    const Matrix tmp3d = args(3).matrix_value();
    min_Rx_d=tmp3d.data()[0]; dx_d=tmp3d.data()[3]; max_Rx_d=tmp3d.data()[6];
    min_Ry_d=tmp3d.data()[1]; dy_d=tmp3d.data()[4]; max_Ry_d=tmp3d.data()[7];
    min_Rz_d=tmp3d.data()[2]; dz_d=tmp3d.data()[5]; max_Rz_d=tmp3d.data()[8];
    Nx = (dream_idx_type) ((max_Rx_d - min_Rx_d)/dx_d+1.0);
    Ny = (dream_idx_type) ((max_Ry_d - min_Ry_d)/dy_d+1.0);
    Nz = (dream_idx_type) ((max_Rz_d - min_Rz_d)/dz_d+1.0);
  }

  dream_idx_type No = Nx*Ny*Nz;

  //
  // Temporal sampling parameter.
  //

  if (!(mxGetM(4) == 1 && mxGetN(4) == 1) ) {
    error("Argument 5 must be a scalar!");
    return oct_retval;
  }

  float dt_f = 0.0;
  double dt_d = 0.0;

  if (use_float) {
    const FloatMatrix tmp4f= args(4).float_matrix_value();
    const float *s_par_f = tmp4f.data();
    dt_f = s_par_f[0];
  } else {
    const Matrix tmp4d = args(4).matrix_value();
    const double *s_par_d = (double*) tmp4d.data();
    dt_d = s_par_d[0];
  }

  //
  // Start point of impulse response vector ([us]).
  //

  if (!ap.check_delay("das_uni", args, 5, No)) {
    return oct_retval;
  }

  DelayType delay_type = DelayType::single;  // delay is a scalar.
  if (mxGetM(5) * mxGetN(5) != 1) {
    delay_type = DelayType::multiple; // delay is a vector.
  }

  const float *delay_f = nullptr;
  const double *delay_d = nullptr;

  if (use_float) {
    const FloatMatrix tmp5f = args(5).float_matrix_value();
    delay_f = tmp5f.data();
  } else {
    const Matrix tmp5d = args(5).matrix_value();
    delay_d = tmp5d.data();
  }


  //
  // Material parameter
  //

  // Check that arg 7 is a scalar.
  if (!(mxGetM(6)==1 && mxGetN(6)==1)) {
    error("Argument 7 must be a scalar!");
    return oct_retval;
  }

  // Sound speed.
  float cp_f = 0.0;
  double cp_d = 0.0;

  if (use_float) {
    const FloatMatrix tmp6f = args(6).float_matrix_value();
    const float *m_par_f = tmp6f.data();
    cp_f = m_par_f[0];
  } else {
    const Matrix tmp6d = args(6).matrix_value();
    const double *m_par_d = (double*) tmp6d.data();
    cp_d = m_par_d[0];
  }

  //
  // DAS method
  //

  DASType das_type;
  if (!ap.parse_das_arg("das_uni", args, 7, das_type)) {
    return oct_retval;
  }

  //
  // Init DAS_UNI and output arg.
  //

  Matrix Im_mat_d;
  FloatMatrix Im_mat_f;
  float *Im_f = nullptr;
  double *Im_d = nullptr;

  bool init_das = true;

  if (use_float) {

    // Create an output matrix for the impulse response.
    Im_mat_f = FloatMatrix(No,1);
    Im_f = (float*) Im_mat_f.data();

    //
    // Single precision DAS_UNI
    //

    if (das_uni_f) { // das_uni (float) object exist - check if we can reuse previous das_uni init
      if (!das_uni_f->das_setup_has_changed(das_type, a_scan_len,
                                            min_t_f, pitch_t_f, max_t_f,
                                            min_r_f, pitch_r_f, max_r_f,
                                            min_Rx_f, dx_f, max_Rx_f,
                                            min_Ry_f, dy_f, max_Ry_f,
                                            min_Rz_f, dz_f, max_Rz_f)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_uni_d) {
        das_uni_d = nullptr; // Release the double object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_uni_f = std::make_unique<DAS_UNI<float>>(das_type, a_scan_len,
                                                     min_t_f, pitch_t_f, max_t_f,
                                                     min_r_f, pitch_r_f, max_r_f,
                                                     min_Rx_f, dx_f, max_Rx_f,
                                                     min_Ry_f, dy_f, max_Ry_f,
                                                     min_Rz_f, dz_f, max_Rz_f);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return oct_retval;
      }
    }

  } else {

    //
    // Double precision DAS_UNI
    //

    // Create an output matrix for the impulse response.
    Im_mat_d = Matrix(No,1);
    Im_d = (double*) Im_mat_d.data();

    if (das_uni_d) { // das_uni (double) object exist - check if we can reuse previous das_uni init
      if (!das_uni_d->das_setup_has_changed(das_type, a_scan_len,
                                            min_t_d, pitch_t_d, max_t_d,
                                            min_r_d, pitch_r_d, max_r_d,
                                            min_Rx_d, dx_d, max_Rx_d,
                                            min_Ry_d, dy_d, max_Ry_d,
                                            min_Rz_d, dz_d, max_Rz_d)) {
        init_das = false;
      }
    }

    if (init_das) {

      if (das_uni_f) {
        das_uni_f = nullptr; // Release the float object if it exist to free, in particular, GPU resources (call destructor).
      }

      try {
        das_uni_d = std::make_unique<DAS_UNI<double>>(das_type, a_scan_len,
                                            min_t_d, pitch_t_d, max_t_d,
                                            min_r_d, pitch_r_d, max_r_d,
                                            min_Rx_d, dx_d, max_Rx_d,
                                            min_Ry_d, dy_d, max_Ry_d,
                                            min_Rz_d, dz_d, max_Rz_d);
      }

      catch (std::runtime_error &err) {
        std::cout << err.what();
        return oct_retval;
      }
    }
  }

  if (use_float) { // Single precision
    das_uni_f->cl_das_uni(Yf,
                          min_t_f, pitch_t_f, max_t_f,
                          min_r_f, pitch_r_f, max_r_f,
                          min_Rx_f, dx_f, max_Rx_f,
                          min_Ry_f, dy_f, max_Ry_f,
                          min_Rz_f, dz_f, max_Rz_f,
                          dt_f,
                          delay_f[0], cp_f,
                          Im_f);
  } else { // Double precision
    das_uni_d->cl_das_uni(Yd,
                          min_t_d, pitch_t_d, max_t_d,
                          min_r_d, pitch_r_d, max_r_d,
                          min_Rx_d, dx_d, max_Rx_d,
                          min_Ry_d, dy_d, max_Ry_d,
                          min_Rz_d, dz_d, max_Rz_d,
                          dt_d,
                          delay_d[0], cp_d,
                          Im_d);
  }

  if (use_float) {
    oct_retval.append(Im_mat_f);
  } else {
    oct_retval.append(Im_mat_d);
  }

  return oct_retval;
}
