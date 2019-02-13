function d = arr_focus_delays(gx,cp,foc,ang)
% d = arr_focus_delays(x,foc,ang)
%
% Function to compute array focusing delays.
%
% x   - Vector of (horizontal) center positions of the array elements.
% cp  - Soundspeed of the media [m/s].
% foc - Focal depth [mm].
% ang - Steering angle [deg].
%
% d  - vector of focusing delays [us].
%
% Fredrik Lingvall 2004-02-17.

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

d = zeros(L,1);
for n=1:L
  d(n) = (r_max - norm([xp(n) zp(n)]' - [gx(n) 0]'))*1000/cp;
end
