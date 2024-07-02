% [H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,err_level,device)
%
% DREAMRECT - Computes the spatial impulse response for a rectangular
% transducer using parallel processing (using threads).
%
%  Observation point(s) ([mm]):
%
%  Ro
%      An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ...  xoN yoN
%      zoN]; where N is the number of observation points.
%
%  Geometrical parameters: geom_par = [a b];
%
%  a
%      x-size of the transducer [mm].
%  b
%      y-size of the transducer [mm].
%
%  Sampling parameters: s_par = [dx dy dt nt];
%
%  dx
%      Spatial x-direction discretization size [mm].
%  dy
%      Spatial y-direction discretization size [mm].
%  dt
%      Temporal discretization period (= 1/sampling freq) [us].
%  nt
%      Length of impulse response vector.
%
%  Start point of SIR:
%
%  delay
%      Scalar delay for all observation points or a vector with
%      individual delays for each observation point [us].
%
%  Material parameters: m_par = [v cp alpha];
%
%  v
%      Normal velocity [m/s].
%  cp
%      Sound velocity [m/s].
%  alpha
%      Attenuation coefficient [dB/(cm MHz)].
%
%  Error Handling: err_level; err_level is an optional text string
%  parameter for controlling the error behavior, options are:
%
%  'ignore'
%      An error is ignored (no error message is printed and the
%      program is not stopped) but the err output argument is
%      non-zero if an error occured.
%  'warn'
%      An error message is printed but the program in not stopped
%      (and err is non-zero).
%  'stop'
%      An error message is printed and the program is stopped.
%  'device'
%      A string which can be one of 'cpu' or 'gpu'.
%
% dreamrect is a mex-function that is a part of the DREAM Toolbox
% available at <https:/ithub.com/frli8848/DREAM>.
%
% Copyright © 2006-2024 Fredrik Lingvall.
