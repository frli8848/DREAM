function [d_vec] = arr_focus_delays(gx,cp,foc,ang)
  %% d = arr_focus_delays(x,foc,ang)
  %%
  %% Function to compute array focusing delays.
  %%
  %% x   - Vector of (horizontal) center positions of the array elements.
  %% cp  - Soundspeed of the media [m/s].
  %% foc - Focal depth [mm].
  %% ang - Steering angle [deg].
  %%
  %% d_vec - vector of focusing delays [us].
  %%
  %% Copyright (C) 2004,2009,2023 Fredrik Lingvall

  L = length(gx);

  xp = foc*tan(ang*pi/180)*ones(L,1);
  zp = foc*ones(L,1);

  r_max = 0;
  for l=1:L
    tmp_max = norm([xp(l) zp(l)]' - [gx(l) 0]');
    if tmp_max > r_max
      r_max = tmp_max;
    end
  end

  d_vec = zeros(L,1);
  for n=1:L
    d_vec(n) = (r_max - norm([xp(n) zp(n)]' - [gx(n) 0]'))*1000/cp;
  end

  %% Virtual source (defocused array)
  if (foc < 0.0)
    mx_foc_us = max(d_vec);
    d_vec = abs(d_vec - mx_foc_us);
  end

end
