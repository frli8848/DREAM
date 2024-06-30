% [Im] = das_uni(Y,gt,gr,Ro,dt,delay,cp,method).
%
% DAS_UNI Computes the delay-and-sum processed reconstruction
% (beamformed image) for three different array geometries: SAFT, TFM,
% and RCA using uniform grids both for the array element positions
% and for the image points.  In the synthtic aperture focusing
% techinque (SAFT) one use one transmitter and one reciever (the
% same) and moves that along the array aperture.  In the total
% focusing technique (TFM) one first transmit with the first element
% and then recieve with all elements, transmit with the second
% element and again receive with all elements and so on until the
% last transmit element (ie., transmit with all - receive with all).
% The last method, row-column adressed (RCA) array, is a variant of
% TFM where the elements are arranged in a crossed layout to form a
% 2D array.
%
% Data matrix:
%
% ‘Y’
%      A K x N matrix where K is the A-scan length and the number of
%      A-scans, N, depends on DAS algorithm selected (see below).
%
% Transmit element vector [mm]:
%
% ‘gt’
%      A 3 element vector, gt = [min_t, pitch_t, max_t].
%
% Receive element vector [mm]:
%
% ‘gr’
%      A 3 element vector, gr = [min_r, pitch_r, max_r].
%
% Observation (image) point matrix [mm]:
%
% ‘Ro’
%      A 3 x 3 matrix, Ro = [min_Rx, dx, max_Rx; min_Ry, dy, max_Ry;
%      min_Rz, dz, max_Rz].
%
% Sampling parameter dt:
%
% ‘dt’
%      Temporal discretization period (= 1/sampling freq) [us].
%
% Data start and pulse delay compensation:
%
% ‘delay’
%      Scalar delay for all observation points or a vector with
%      individual delays for each observation point [us].
%
% Sound speed:
%
% ‘cp’
%      Sound velocity of the medium [m/s].
%
% DAS algorithm: das_met is a text string parameter for selecting DAS
% algorithm, options are:
%
% ‘'saft'’
%      When SAFT is selected data Y must be an K x Lt (SAFT is also
%      selected if Gr=[]).
% ‘'tfm'’
%      When TFM is selected a linear array is assumed and data Y must
%      be a K x Lt*Lr matrix.
% ‘'rca_coltx' or 'rca_coltx'’
%      When RCA is selected a 2D RCA array is assumed and data Y must
%      be a K x Lt*Lr matrix.
%
% das_uni is a mex-function that is a part of the DREAM Toolbox
% available at <https://github.com/frli8848/DREAM>.
%
% Copyright © 2023 Fredrik Lingvall.
