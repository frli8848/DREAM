function [cp] = h2o_soundspeed(T,Pin,unit,S)
%% [cp] = h2o_soundspeed(T,P,unit,S)
%%
%% Computes the sound speed of (sea) water as a function of
%% temperature, salinity, and pressure. The sound speed is
%% computed using a formula from:
%%
%% Del Grosso, V. A., "New equation for the speed of sound in natural
%%    waters (with comparisons to other equations)", J. Acoust. Soc. Am., Vol.56,
%%    No. 4, October 1974, pp. 1084--1091.
%%
%% Input parameters:
%%
%% T    - Temperture [in degrees Celsius].
%% P    - Pressure (optional, default 1 [atm]).
%% unit - Text string defining the pressure unit of arg 2
%%        ['Pa','bar','at','atm','mmHg', or, 'psi']  (optional,
%%        default ['Pa']).
%% S    - Salinity [in parts per thousand (ppt)] (optional, default 0 [ppt]).
%%
%% Copyright (C) 2007,2009,2023 Fredrik Lingvall

  if nargin < 4
    S = 0;
  end

  if nargin < 3
    unit = 'Pa';                % Default to the SI unit.
  end

  if nargin < 2
    Pin = 101325;               % Default to 1 [atm].
  end

  %% Conversion table from: http://en.wikipedia.org/wiki/Pressure

  switch lower(unit)

    case {'pa'}                 % Pascal.
      P = Pin * 10.197e-6;

    case {'bar'}                % Bar.
      P = Pin * 1.0197;

    case {'at'}                 % Technical atmosphere.
      P = Pin;

    case {'atm'}                % Atmosphere.
      P = Pin * 1.0332;

    case {'mmhg'}               % Torr.
      P = Pin * 1.3595e-3;

    case {'psi'}                % Pound-force per square inch.
      P = Pin * 70.307e-3;

    otherwise
      disp('Unknown pressure unit - using Pascal!')
      P = Pin * 10.197e-6;
end

  if nargin < 1,
    help h2o_soundspeed;
    return;
  end

  C_000 = 1402.392;

  deltaC_T = 0.501109398873e1 .* T ...
             - 0.550946843172e-1 .* T.^2 ...
             + 0.221535969240e-3  .* T.^3;

  deltaC_S = 0.132952290781e1   .* S ...
             + 0.128955756844e-3 .* S.^3;

  deltaC_P = 0.156059257041     .* P ...
             + 0.244998688441e-4 .* P.^2 ...
             - 0.883392332513e-8 .* P.^3;

  deltaC_STP = - 0.127562783426e-1.*T.*S  ...
               + 0.635191613389e-2 .* T    .* P    ...
               + 0.265484716608e-7 .* T.^2 .* P.^2 ...
               - 0.159349479045e-5 .* T    .* P.^3 ...
               + 0.522116437235e-9 .* T    .* P.^3 ...
               - 0.438031096213e-6 .* T.^3 .* P    ...
               - 0.161674495909-8  .* S.^2 .* P.^2 ...
               + 0.968403156410e-4 .* T.^2 .* S    ...
               + 0.485639620015e-5 .* T    .* S.^2 .* P ...
               - 0.340597039004e-3 .* T    .* S    .* P;

  cp = C_000 + deltaC_T + deltaC_S + deltaC_P + deltaC_STP;
end
