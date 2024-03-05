%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Make a plot of the speed for the convolution functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN = [];
T0 = [];
T00 = [];
T1 = [];
T2 = [];
T3 = [];

% Average M times.
M = 10;

disp('This may take some time - go and make some tea...');

fprintf(['\n---------------------------------------------------------' ...
         '-----------------------|\n']);
fprintf(' n\t| conv\t\t|  conv_p\t| fftconv_p\t| fftconv_p wisd.\t|\n');
fprintf(['---------------------------------------------------------' ...
         '-----------------------|\n']);
for n=[1 10:10:100 200:100:2000]

  L1 = n;
  L2 = n;
  K=100;
  X = randn(L1,K);
  Y = randn(L2,K);

  % Reference conv method.
  t0 = 0;
  for m=1:M
    eval('tic');
    for k=1:K
      Z0 = conv(X(:,k),Y(:,k));
    end
    t0 = t0 + toc;
  end
  t0 = t0/M;

  % The fftconv function is only available
  % in Octave.
  if exist('OCTAVE_VERSION')
    t00 = 0;
    for m=1:M
      eval('tic');
      for k=1:K
        Z00 = fftconv(X(:,k),Y(:,k));
        %Z00 = conv(X(:,k),Y(:,k));
      end
      t00 = t00 + toc;
    end
    t00 = t00/M;
  end

  t1 = 0;
  for m=1:M
    eval('tic');
    Z = conv_p(X, Y);
    t1 = t1 + toc;
  end
  t = t1/M;

  %clear wisdom
  t2 = 0;
  for m=1:M
    eval('tic');
    [Z2,wisdom] = fftconv_p(X, Y);
    t2 = t2+toc;
  end
  t2 = t2/M;

  t3 = 0;
  for m=1:M
    eval('tic');
    Z3 = fftconv_p(X, Y, wisdom);
    t3 = t3 + toc;
  end
  t3 = t3/M;

  fprintf('%d\t| %f\t| %f\t| %f\t| %f\t\t|\n',n,t0,t1,t2,t3);

  if norm(Z-Z2) > 1e-10
    disp('Warning: large error in conv_p/fftconv_p!');
  end

  NN =  [NN n];
  T0 = [T0 t0];
  if exist('OCTAVE_VERSION')
    T00 = [T00 t00];
  end
  T1 = [T1 t1];
  T2 = [T2 t2];
  T3 = [T3 t3];
end

figure(1);
clf

semilogy(NN,T1)
hold('on');
semilogy(NN,T2,'r');
semilogy(NN,T3,'g');
semilogy(NN,T0,'.-');
if exist('OCTAVE_VERSION')
  semilogy(NN,T00,':');
end

hold('off');
xlabel('Vector length');
ylabel('Time [s]');
grid('on');

if exist('OCTAVE_VERSION')
  legend('conv\_p','fftconv\_p', 'fftconv\_p w. wisd.','conv','fftconv');
else
  legend('conv\_p','fftconv\_p', 'fftconv\_p w. wisd.','conv');
end
