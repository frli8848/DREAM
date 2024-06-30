% DREAM Toolbox functions
%
% DR-based transducer functions:
%
%    dreamline         - Strip/line transducer
%    dreamrect         - Rectangular transducer
%    dreamrect_f       - Focused rectangular transducer
%    dreamcirc         - Circular transducer
%    dreamcirc_f       - Focused circular transducer
%    dreamsphere       - Spherical concave/convex transducer (defocused/defocused)
%    dreamcylind       - Cylindrical concave/convex transducer (focused/defocused)
%    dream_arr_rect    - Array with rectangular elements
%    dream_arr_circ    - Array with circular elements
%    dream_arr_cylind  - Array with cylindrical concave/convex elements (focused/defocused)
%    dream_arr_annu    - Annular array
%
% Analytic transducer functions:
%
%    rect_sir  -  time-continous (analytic) rectangual ransducer
%    circ_sir  - time-continous (analytic) circular ransducer
%    scirc_sir - sampled time-continous (analytic) circular transducer
%
% Threaded (parallel) signal processing algorithms:
%
%    conv_p      - threaded one dimensional convolution
%    fftconv_p   - threaded one dimensional FFT based convolutions
%    fftconv_ola - threaded one dimensional FFT based convolutions using overlap-and-add
%    sum_fftconv - threaded sum of one dimensional FFT based convolutions
%    das         - threaded delay-and-sum (DAS) beamforming (SAFT, TFM, and RCA)
%    das_uni     - (OpenCL only) DAS beamforming using uniform grids
%
% Copyright (C) 2006-2024 Fredrik Lingvall.
