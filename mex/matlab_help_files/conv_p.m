% [Y] = conv_p(A, B).
%
% CONV_P Computes one dimensional convolutions of the columns in the
% matrix A and the matrix (or vector) B.
%
% Input parameters:
%
% ‘A’
%  An MxN matrix.
% ‘B’
%  A KxN matrix or a K-length vector.  If B is a vector each
%  column in A is convolved with the vector B.
%
%  Output parameter:
%
%  ‘Y’
%   The (M+K-1)xN output matrix.
%
% conv_p is a mex-function that is a part of the DREAM Toolbox
% available at <https:/ithub.com/frli8848/DREAM>.
%
% Copyright © 2006-2023 Fredrik Lingvall.
