function plot_das_fmc_delays(Y,Gt,Gr,Ro,dt,delay,cp,n)
  %% plot_das_fmc_delays(Y,Gt,Gr,Ro,dt,delay,cp,n)
  %%
  %% A function for investigating the delay-and-sum TFM
  %% type (linear array and RCA array) the delays
  %% and the pulse compensations delay.
  %%
  %% It can be used with the test_das_tfm test_das_rca
  %% scripts like for example,
  %%
  %% >> DO_PLOTTING=1;
  %% >> test_das_tfm
  %% >> clf
  %% >> a_scan_idx = 1;
  %% >> plot_das_fmc_delays(Yfmc,Gt,Gr,[0 0 10],dt,delay,cp,a_scam_idx);
  %%
  %% where we plot the 1st A-scan in the simulated FMC data, Yfmc,
  %% and the corresponding TFM and FMC focusing delays.

  xo = Ro(1);
  yo = Ro(2);
  zo = Ro(3);

  KN = size(Y,2);
  assert(KN == size(Gt,1)*size(Gr,1));

  [ind_t,ind_r] =  ind2sub([size(Gt,1) size(Gr,1)], n)

  delay_ms = delay/1000.0;

  %%
  %% TFM (linear array) DAS delay
  %%

  %% Transmit
  gx_t = Gt(ind_t,1) - xo;
  gy_t = Gt(ind_t,2) - yo;
  gz_t = Gt(ind_t,3) - zo;

  %% Receive
  gx_r = Gr(ind_r,1) - xo;
  gy_r = Gr(ind_r,2) - yo;
  gz_r = Gr(ind_r,3)- zo;

  t_tfm = sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t)/cp + ...
          sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r)/cp + ...
          delay_ms;

  t_idx_tfm = t_tfm/dt*1000.0;

  [gx_t gy_t  gz_t; gx_r gy_r gz_r]

  %%
  %% TFM RCA DAS delay
  %%

  %% Element (stripe) lengths
  gt_y_min = Gr(1,2);
  gt_y_max = Gr(end,2);
  gr_x_min = Gt(1,1);
  gr_x_max = Gt(end,1);

  %% Transmit
  gx_t = Gt(ind_t,1) - xo;
  if ( (yo >= gt_y_min) && yo <= gt_y_max)
    gy_t = 0.0;               % We are inside the stripe aperture.
  else
    if (yo < gt_y_min)
      gy_t = gt_y_min - yo;
    else
      gy_t = gt_y_max - yo;
    end
  end
  gz_t = Gt(ind_t,3)- zo;

  %% Receive
  gx_r = Gr(ind_r,1) - xo;
  if ( (yo >= gr_x_min) && yo <= gr_x_max)
    gy_r = 0.0;               % We are inside the stripe aperture.
  else
    if (yo < gr_x_min)
      gy_r = gr_x_min - yo;
    else
      gy_r = gr_x_max - yo;
    end
  end
  gy_r = Gr(ind_r,2) - yo;
  gz_r = Gr(ind_r,3) - zo;

  t_rca = sqrt(gx_t*gx_t + gy_t*gy_t + gz_t*gz_t)/cp + ...
          sqrt(gx_r*gx_r + gy_r*gy_r + gz_r*gz_r)/cp + ...
          delay_ms;

  t_idx_rca = t_rca/dt*1000.0;

  [gx_t gy_t  gz_t; gx_r gy_r gz_r]

  plot(Y(:,n));
  ax = axis;
  hold on
  l1 = line([t_idx_tfm  t_idx_tfm],[ax(3) ax(4)]);
  set(l1,'color','red')
  l2 = line([t_idx_rca  t_idx_rca],[ax(3) ax(4)]);
  set(l2,'color','green')
end
