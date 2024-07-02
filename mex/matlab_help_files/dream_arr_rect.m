% [H,err] = dream_arr_rect(Ro,geom_par,G,s_par,delay,m_par,foc_met,...
%     focal,steer_met,steer_par,apod_met,apod,win_par,err_level);
%
%     DREAM_ARR_RECT - Computes the spatial impulse response for an array
%     transducer with rectangular elements using parallel processing
%     (using threads).
%
%     Observation point(s) ([mm]):
%
%     Ro
%          An N x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ...  xoN yoN
%          zoN]; where N is the number of observation points.
%
%     Geometrical parameters: geom_par = [a b];
%
%     a
%          Element width in x-direction [mm].
%     b
%          Element width in y-direction [mm].
%
%     Array grid parameter:
%
%     G
%          An L x 3 grid function matrix.  Column 1 contain the
%          x-postions of the elements, column 2 the y-positions, and
%          column 3 the z-positions, and where L is the number of
%          elements.
%
%     Sampling parameters: s_par = [dx dy dt nt];
%
%     dx
%          Spatial x-direction discretization size [mm].
%     dy
%          Spatial y-direction discretization size [mm].
%     dt
%          Temporal discretization period (= 1/sampling freq) [us].
%     nt
%          Length of impulse response vector.
%
%     Start point of SIR:
%
%     delay
%          Scalar delay for all observation points or a vector with
%          individual delays for each observation point [us].
%
%     Material parameters: m_par = [v cp alpha];
%
%     v
%          Normal velocity [m/s].
%     cp
%          Sound velocity [m/s].
%     alpha
%          Attenuation coefficient [dB/(cm MHz)].
%
%     Focusing parameters: foc_met and focal:
%
%     foc_met
%          Focusing method, options are: 'off', 'x', 'y', 'xy', 'x+y',
%          and 'ud'.
%     focal
%          Focal distance [mm].  If foc_met = 'ud' (user defined) then
%          focal is a vector of focusing delays.
%
%     Beam steering parameters: steer_met and steer_par:
%
%     steer_met
%          Beam steering method, options are: 'off', 'x', 'y', and 'xy'.
%     steer_par = [theta phi];
%          theta [deg] is the x-direction steer angle and phi [deg] is
%          the y-direction steer angle.
%
%     Apodization parameters: apod_met, apod, and win_par.  The apod_met
%     (apodization method) options are:
%
%     'off'
%          No apodization.
%     'ud'
%          User defined apodization.
%     'triangle'
%          Triangle window.
%     'gauss'
%          Gaussian (bell-shaped) window.
%     'raised'
%          Raised cosine.
%     'hann'
%          Hann window.
%     'hamming'
%          Hamming window.
%     'simply'
%          Simply supported.
%     'clamped'
%          Clamped.
%
%     and the apod and win_par parameters are:
%
%     apod
%          Vector of apodiztion weights (used for the 'ud' option).
%     win_par
%          A scalar parameter for raised cosine and Gaussian apodization
%          functions.
%
%     Error Handling: err_level; err_level is an optional text string
%     parameter for controlling the error behavior, options are:
%
%     'ignore'
%          An error is ignored (no error message is printed and the
%          program is not stopped) but the err output argument is
%          non-zero if an error occured.
%     'warn'
%          A warning message is printed but the program in not stopped.
%     'stop'
%          An error message is printed and the program is stopped.
%
% dream_arr_rect is a mex-function that is a part of the DREAM
% Toolbox available at <https:/ithub.com/frli8848/DREAM>.
%
% Copyright © 2006-2024 Fredrik Lingvall.
