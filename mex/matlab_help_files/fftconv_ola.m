% [Y,wisdom_str_out] = fftconv_ola(A,B,block_len,wisdom_str_in);
%
%     FFTCONV_OLA - Computes the one dimensional convolution of the
%     columns in matrix A and the matrix (or vector) B using the
%     overlap-and-add method.
%
%     Input parameters:
%
%     ‘A’
%          An M x N matrix.
%     ‘B’
%          A K x N matrix or a K-length vector.  If B is a vector each
%          column in A is convolved with the vector B.
%     ‘wisdom_str_in’
%          Optional parameter.  If the wisdom_str_in parameter is not
%          supplied then ‘fftconv_ola’ calls fftw wisdom plan functions
%          before performing any frequency domain operations.  This
%          overhead can be avoided by supplying a pre-computed fftw
%          wisdom string wisdom_str_in.  For more information see the
%          fftw user manunal available at <http://www.fftw.org>.
%
%     Output parameters:
%
%     ‘Y’
%          The (M+K-1) x N output matrix.
%     ‘wisdom_str_out’
%          Optional parameter.  If the wisdom_str_out output parameter is
%          supplied then ‘fftconv_ola’ will call fftw wisdom plan
%          functions and return the wisdom string which then can be used
%          to speed up subsequent calls to ‘fftconv_ola’ by suppying the
%          string as the input argument wisdom_str_in.
%
%     In place modes:
%
%     fftconv_ola(A,B,block_len,Y);
%
%     fftconv_ola(A,B,block_len,Y,mode);
%
%     fftconv_ola(A,B,block_len,Y,wisdom_str_in);
%
%     fftconv_ola(A,B,block_len,Y,mode,wisdom_str_in);
%
%     where the 'mode' is a string which can be '=', '+=', or '-='.
%
%     NOTE: fftconv_ola requires the FFTW library version 3
%     <http://www.fftw.org>.
%
% fftconv_ola is a mex-function that is a part of the DREAM Toolbox
% available at <https:/ithub.com/frli8848/DREAM>.
%
% Copyright © 2010-2021 Fredrik Lingvall.
