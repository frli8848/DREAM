% [Im] = das(Y,Gt,Gr,Ro,dt,delay,cp,method,err_level,device).
%
% DAS Computes the delay-and-sum processed reconstruction (beamformed
% image) for three different array geometries: SAFT, TFM, and RCA. In
% the synthtic aperture focusing techinque (SAFT) one use one
% transmitter and one reciever (the same) and moves that along the
% array aperture.  In the total focusing technique (TFM) one first
% transmit with the first element and then recieve with all elements,
% transmit with the second element and again receive with all
% elements and so on until the last transmit element (ie., transmit
% with all - receive with all).  The last method, row-column adressed
% (RCA) array, is a variant of TFM where the elements are arranged in
% a crossed layout to form a 2D array.
%
% Data matrix:
%
% Y
%      A K x N matrix where K is the A-scan length and the number of
%      A-scans, N, depends on DAS algorithm selected (see below).
%
% Transmit element grid matrix:
%
% Gt
%      An Lt x 3 matrix, Gt = [x1 y1 z2; x2 y2 z2; ...  xLt yoLt
%      zLt]; where Lt is the number of tranmsit elements.
%
% Receive element grid matrix:
%
% Gr
%      An Lr x 3 matrix, Gr = Gt = [x1 y1 z2; x2 y2 z2; ...  xLr yoLr
%      zLr]; where Lr is the number of tranmsit elements.
%
% Observation point(s) ([mm]):
%
% Ro
%      An No x 3 matrix, Ro = [xo1 yo1 zo2; xo2 yo2 zo2; ...  xoNo
%      yoNo zoNo]; where No is the number of observation points.
%
% Sampling parameter dt:
%
% dt
%      Temporal discretization period (= 1/sampling freq) [us].
%
% Data start and pulse delay compensation:
%
% delay
%      Scalar delay for all observation points or a vector with
%      individual delays for each observation point [us].
%
% Sound speed:
%
% cp
%      Sound velocity of the medium [m/s].
%
% DAS algorithm: das_met is a text string parameter for selecting DAS
% algorithm, options are:
%
% 'saft'
%      When SAFT is selected data Y must be an K x Lt (SAFT is also
%      selected if Gr=[]).
% 'tfm'
%      When TFM is selected a linear array is assumed and data Y must
%      be a K x Lt*Lr matrix.
% 'rca_coltx' or 'rca_rowtx'
%      When RCA is selected a 2D RCA array is assumed and data Y must
%      be a K x Lt*Lr matrix.
%
% Error Handling: err_level is an optional text string parameter for
% controlling the error behavior, options are:
%
% 'ignore'
%      An error is ignored (no error message is printed and the
%      program is not stopped) but the err output argument is
%      non-zero if an error occured.
% 'warn'
%      An error message is printed but the program in not stopped
%      (and err is non-zero).
% 'stop'
%      An error message is printed and the program is stopped.
%
% Compute device:
%
% 'device'
%      A string which can be one of 'cpu' or 'gpu'.
% 'verbose'
%      A string which when it is 'verbose' will printout OpenCL device info.
%
% das is a mex-function that is a part of the DREAM Toolbox
% available at <https://github.com/frli8848/DREAM>.
%
% Copyright © 2008-2025 Fredrik Lingvall.
