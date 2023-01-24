%%
%% ------------- Rectangular Transducer --------------------------
%%

setup_dream_parameters;

% Geometrical parameters.
a = 10;				% x-size.
b = 15;				% y-size.
geom_par = [a b];

%% Make sure the number of obs points are a 64 (the OpenCL workgroup size)
[X,Y,Z] = meshgrid(-15:0.5:15,-15:0.5:15,1.25:0.25:65);
Ro =[X(:) Y(:) Z(:)];
fprintf('Num obs points: %d\n', size(Ro,1))

tic,[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop');toc

tic,[H,err] = dreamrect(Ro,geom_par,s_par,delay,m_par,'stop','gpu');toc
