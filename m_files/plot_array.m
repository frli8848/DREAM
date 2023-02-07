function plot_array(gx,a,foc,cp)
  %% plot_array(gx,a,foc,cp)
  %%
  %% The observation points do not need to be regurlary sampled or
  %% ordered in any particular way.
  %%
  %% Input arguments:
  %% gx       - Position vector for the array elements [mm].
  %% a        - element size [mm].
  %% foc      - Focal depth or focusing delay vector.
  %% theta    - Steering angle (degrees).
  %%
  %% Copyright (C) 2009,2023 Fredrik Lingvall.

  %% Number of array elements.
  Lt = size(gx,1);
  Lr = size(gx,1);

  %% Make the array casing 5 mm high.
  y_max = 0;
  y_min = -5;

  %% Draw outline of array.
  h = line([min(gx)-a/2 max(gx)+a/2], [y_min y_min]);
  set(h,'LineWidth',2);
  h = line([min(gx)-a/2 max(gx)+a/2],[y_max y_max]);

  h = line([min(gx)-a/2  min(gx)-a/2],[y_min y_max]);
  set(h,'LineWidth',2);
  h = line([max(gx)+a/2 max(gx)+a/2],[y_min y_max]);
  set(h,'LineWidth',2);

  for l=1:Lr
    el_pos = gx(l);

    line([el_pos-a/2 el_pos+a/2], [y_min y_min]);
    line([el_pos-a/2 el_pos+a/2 ],[y_max y_max]);

    line([el_pos-a/2 el_pos-a/2],[y_min y_max]);
    line([el_pos+a/2 el_pos+a/2],[y_min y_max]);
  end

  if Lt>1
    for l=1:Lt
      el_pos = gx(l);

      y_m = y_max - foc(l)*cp/1e3;

      fill([el_pos-a/2 el_pos-a/2 el_pos+a/2 el_pos+a/2], ...
           [y_min y_m  y_m y_min], ...
           [0.5 0.5 0.5]);
    end

  else
    el_pos = trans(1,1);

    fill([el_pos-a/2 el_pos-a/2 el_pos+a/2 el_pos+a/2], ...
         [y_min y_max  y_max y_min], ...
         [0.5 0.5 0.5]);
end
