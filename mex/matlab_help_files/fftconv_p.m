% [Y,wisdom_str_out] = fftconv_p(A,B,wisdom_str_in);
%
%    FFTCONV_P - Computes the one dimensional convolution of the columns
%    in matrix A and the matrix (or vector) B.
%
%    Input parameters:
%
%    ‘A’
%         An M x N matrix.
%    ‘B’
%         A K x N matrix or a K-length vector.  If B is a vector each
%         column in A is convolved with the vector B.
%
%    ‘wisdom_str_in’
%
%         Optional parameter.  If the wisdom_str_in parameter is not
%         supplied then ‘fftconv_p’ calls FFTW wisdom plan functions
%         before performing any frequency domain operations.  This
%         overhead can be avoided by supplying a pre-computed FFTW
%         wisdom string wisdom_str_in.  For more information see the
%         FFTW user manunal available at <http://www.fftw.org>.
%
%    Output parameters:
%
%    ‘Y’
%         The (M+K-1) x N output matrix.
%
%    ‘wisdom_str_out’
%         Optional parameter.  If the wisdom_str_out output parameter is
%         supplied then ‘fftconv_p’ will call FFTW wisdom plan functions
%         and return the wisdom string which then can be used to speed
%         up subsequent calls to ‘fftconv_p’ by suppying the string as
%         the input argument wisdom_str_in.
%
%         In place modes:
%
%         fftconv_p(A,B,Y);
%         fftconv_p(A,B,Y,mode);
%         fftconv_p(A,B,Y,wisdom_str_in);
%         fftconv_p(A,B,Y,mode,wisdom_str_in);
%
%         where the 'mode' is a string which can be '=', '+=', or '-='.
%
% fftconv_p is a part of the DREAM Toolbox available at
% <https:/ithub.com/frli8848/DREAM>.
%
% Copyright © 2006-2023 Fredrik Lingvall.
