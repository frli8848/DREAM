%
% Fredrik Lingvall 2002-11-08
% Fredrik Lingvall 2003-09-15
% Fredrik Lingvall 2005-10-18
% Fredrik Lingvall 2007-02-12


%
% Front page.
%

PRINT=0;

if 1
  
  if 1
    
    Fs = 10;
    Fs = Fs*1e6;
    Ts = 1/Fs;
    c_H2O = 1484;
    %c_H2O = 1000;
    
    %z = 190;
    z = 10;
    t_z = z*1e-3/c_H2O * 1e6; % us 
    
    % Observation point.
    xo = 0;
    yo = 0;
    zo = z;
    ro = [xo yo zo];
    
    % Descretization parameters.
    dx = 0.01; % mm.
    dy = 0.01; % m.m
    dt = Ts*1e6; % us. 
    
    
    % Length of spatial impulse response vector.
    %nt = 1000;
    nt = 600;
    
    t = 0:Ts:Ts*(nt-1);
    
    % Material parameters.
    v     = 1.0; 				% Normal velocity.
    cp    = c_H2O; 				% Sound speed.
    alfa  = 0; 				% Absorbtion (dB/cm Hz).
    %dens  = 1000; 				% Density.
    
    % Grid function (position vectors of the elements).
    %x = -16:1:16;
    %[gx,gy] = meshgrid(x);
    %gx = gx(:);
    %gy = gy(:);
    %gz = zeros(length(gx),1);
    
    gx = 0;
    gy = 0;
    gz = zeros(length(gx),1);
    
    
    
    % ------------- Focused Cylindrical Transducer --------------------------
    
    % Geometrical parameters.
    a = 20;				% x-width of the transducer.
    b = 20;				% y-width of the transducer.
    r = 100;				% Radius of the curvature.
    geom_par = [a b r];
    
    % Sampling parameters.
    s_par = [dx dy dt nt];
    
    % Start point of SIR
    %delay = t_z;
    delay = 0;
    
    % Material parameters.
    m_par = [v cp alfa];
    
    n = 1;
    H = [];  
    %for d=-50:0.5:50
    % xo = d;
    %h = dreamcylind_f([xo yo zo],geom_par,s_par,delay,m_par,'stop');
    %H(:,n)  = h(:);
    % n = n+1;
    %end
    
    xo = -50:0.5:50;
    yo = zeros(size(xo));
    zo = z*ones(size(xo));
    Ro = [xo(:) yo(:) zo(:)];
    H = dreamcylind_f(Ro,geom_par,s_par,delay,m_par,'ignore');
    %H = dreamrect(Ro,[a b],s_par,delay,m_par,'ignore');
    
    [I,J] = find(H == 0);
    
    figure(1)
    clf
    %imagesc(H);
    mesh(H)
    axis tight
    ax = axis;
    
    H2 = H;
    for n=1:length(I)
      H2(I(n),J(n)) = NaN;
    end
    X=1:size(H,2);
    Y=1:size(H,1);
    
    %sp = mesh(H)
    %sp =surfl(X,Y,H,[0 1 0],[0 1 0 0])
    %sp =surfl(H)
    sp =surf(H2)
    %colormap gray
    %colormap hot
    colormap jet
    %colormap summer
    %colormap winter
    axis(ax);
    axis off
    view([135 32])
    %set(sp,'FaceColor','grey');
    %set(sp,'EdgeColor','none');
    set(sp,'FaceLighting','phong');
    light('Position',[0 1 0],'Style','infinite');
    shading interp;
    
    if PRINT
      %print -depsc frontpage.eps
      print -depsc banner.eps
    end
    drawnow
    
    figure(2)
    clf

    subplot 221
    plot(H(:,round(size(H,2)/2)));
    %title('Focused cylindrical transducer')
    %xlabel('t [\mus]')
    axis off
    
    subplot 222
    plot(H(:,round(size(H,2)/3)));
    axis off
    
    subplot 223
    plot(H(:,round(size(H,2)/5)));
    axis off
    
    subplot 224
    plot(H(:,round(size(H,2)/10)));
    axis off
	
    if PRINT
      print -depsc banner2.eps
    end
  end
end

%
% Apodization windows
%

if 1

  n=-0.5:1/1000:0.5;
  N = length(n);
  
  % XTick vector 
  xt = -0.5:0.1:0.5;
  
  figure(4)
  clf
  w = dream_apodwin('triangle',N,1);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  axis([-0.5 0.5 0 1]);
  grid
  
  xlabel('Normalized width','FontSize',18)
  
  if PRINT
    print -depsc eps/triangle.eps
  end
  
  figure(5)
  clf
  
  w = dream_apodwin('gauss',N,0.1);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.3,0.92,'p=0.1');
  set(h,'FontSize',16);
  
  hold

  w = dream_apodwin('gauss',N,1);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.3,0.75,'p=1');
  set(h,'FontSize',16);
  
  w = dream_apodwin('gauss',N,10);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.16,0.5,'p=10');
  set(h,'FontSize',16);
  
  w = dream_apodwin('gauss',N,100);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.07,0.2,'p=100');
  set(h,'FontSize',16);
  
  axis([-0.5 0.5 0 1]);
  grid
  
  xlabel('Normalized width','FontSize',18)
  
  if PRINT
    print -depsc eps/gauss.eps
  end
  
  figure(3)
  clf
  
  w = dream_apodwin('raised',N,0);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.35,-0.5,'p=0');
  set(h,'FontSize',16);
  
  hold

  w = dream_apodwin('raised',N,0.5);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.35,0,'p=0.5')
  set(h,'FontSize',16);
  
  w = dream_apodwin('raised',N,1);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  h=text(0.35,0.5,'p=1');
  set(h,'FontSize',16);
  
  
  axis([-0.5 0.5 -1 2]);
  grid
  xlabel('Normalized width','FontSize',18)
  
  print -depsc eps/raised.eps

  figure(6)
  clf
  
  w = dream_apodwin('simply',N,0);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  axis([-0.5 0.5 0 1]);
  grid

  xlabel('Normalized width','FontSize',18)
  
  if PRINT
    print -depsc eps/simply.eps
  end
  
  figure(7)
  clf
  
  w = dream_apodwin('clamped',N,0);  
  plot(n,w);
  set(gca,'FontSize',16);
  set(gca,'XTick',xt);
  axis([-0.5 0.5 0 1]);
  grid
  xlabel('Normalized width','FontSize',18)

  if PRINT
    print -depsc eps/clamped.eps
  end
  
end
