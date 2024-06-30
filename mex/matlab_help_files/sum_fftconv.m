% Y = sum_fftconv(H, U, wisdom_str);
%
%      SUM_FFTCONV - Computes (using parallel threaded processing) the sum
%      of one dimensional convolutions of the columns in each 2D matrix in
%      the 3D matrix H with the corresponding columns in the matrix U.
%
%      In normal mode sum_fftconv performs an operation similar to:
%
%      fft_len = M+K-1;
%      YF = zeros(fft_len,N);
%      for l=1:L
%        for n=1:N
%          YF(:,n) = YF(:,n) + fft(H(:,n,l),fft_len).* fft(U(:,l),fft_len);
%        end
%      end
%      Y = real(ifft(YF))/fft_len;
%
%      using threaded processing.  The computations are performed using
%      FFT:s.
%
%      Input parameters:
%
%      ‘H’
%           An MxNxL 3D matrix.
%      ‘U’
%           H KxL matrix.
%      ‘wisdom_str’
%           Optional parameter.  If the wisdom_str parameter is not
%           supplied then fftconv calls fftw wisdom plan functions before
%           performing any frequency domain operations.  This overhead can
%           be avoided by supplying a pre-computed fftw wisdom string.
%           For more information see the fftw user manunal available at
%           <http://www.fftw.org>.
%
%      The wisdom_str can be obtained using the fftconv_p function.  A
%      typical example is,
%
%       % Compute a new fftw wisdom string.
%      [tmp,wisdom_str] = fftconv_p(H(:,1,1),U(:,1));
%
%      for i=1:N
%
%        % Do some stuff here.
%
%        Y = sum_fftconv(H, U, wisdom_str);
%      end
%
%      where the overhead of calling fftw plan functions is now avoided
%      inside the for loop.
%
%      Output parameter:
%
%      ‘Y’
%           The (M+K-1)xN output matrix.
%
% In-place mode
% =============
%
% In in-place mode sum_fftconv performs the operations in-place on a
% pre-allocated matrix:
%
%  sum_fftconv(H, U, Y, wisdom_str);
%
% Here sum_fftconv do not have any output arguments and the results are
% instead stored directly in the pre-allocated (M+K-1)xN input matrix Y. A
% typical usage is:
%
% Y = zeros(M+K-1,N); % Allocate space for Y.
%
% for i=1:N
%
%   % Do some stuff here.
%
%    sum_fftconv(H, U, Y, wisdom_str);
% end
%
% where memory allocation of the (possible large) matrix Y now is avoided
% inside the for loop.
%
%    NOTE: A side-effect of using in-place mode is that if a copy Y2 of Y
% is made then both Y2 and Y will be altered by sum_fftconv.  That is, by
% performing,
%
% Y  = zeros(M+K-1,N);
% Y2 = Y; % No actual copy of data here.
% sum_fftconv(H,U,Y,wisdom_str);
%
% then both Y and Y2 will be changed (since Octave do not make a new copy
% of the data in Y2 unless Y2 is changed before the sum_fftconv call).
%
%    sum_fftconv is a part of the DREAM Toolbox available at
%    <https://github.com/frli8848/DREAM>.
%
%    Copyright © 2006-2023 Fredrik Lingvall.
