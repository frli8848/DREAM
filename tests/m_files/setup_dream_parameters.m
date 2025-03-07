ONE_POINT = 1;
%%ONE_POINT = 0; %Many points.
%% DO_PLOTTING = 1;

Fs = 4;   % Sampling freq. [MHz].
Ts = 1/Fs; % [us].

%%
%% Observation point(s).
%%

z = 10;

if ONE_POINT
  % One point.
  xo = 0;
  yo = 0;
  zo = z;
  ro = [xo yo zo];
  Ro = ro;
  %%disp(['Observation point (x,y,z) = ' num2str(Ro)]);
else
  %  Points along x-axis.
  d  = 1;
  xo = (0:d:50);
  yo = zeros(length(xo),1);
  zo = z*ones(length(xo),1);
  Ro = [xo(:) yo(:) zo(:)];
end

%% Descretization parameters.
dx = 0.05;                % [mm].
dy = 0.05;                % [mm]
dt = Ts;                  % [us].
nt = 300;                 % Length of spatial impulse response vector.
s_par = [dx dy dt nt];

t = 0:Ts:Ts*(nt-1);

% Material parameters.
v     = 1.0;                    % Normal velocity.
cp    = 1000;                   % Sound speed.
alpha  = 0.0;                   % Absorbtion (dB/cm Hz).
m_par = [v cp alpha];

t_z = z*1e-3/cp * 1e6;          % us

delay = 0;
%%delay = t_z;
