%
% Script to make the DR-figure for the DREAM user manual.
%

% Fredrik Lingvall 2007-03-20

figure(1);
clf;

s=0:0.01:1;
h_c = plot3(cos(2*pi*s),sin(2*pi*s),zeros(length(s),1));
set(h_c,'LineWidth',2);

hold on;

%h_square = plot3([-0.3 -0.2 -0.2 -0.2 -0.2 -0.3 -0.3 -0.3], [0.3 0.3 0.2 0.2 ...
%                    0.2 0.2 0.2 0.3],zeros(8,1)); 

h_square = patch([-0.3 -0.2 -0.2 -0.2 -0.2 -0.3 -0.3 -0.3], [0.3 0.3 0.2 0.2 ...
                    0.2 0.2 0.2 0.3],zeros(8,1),[0.5 0.5 0.5]);

h_x_arrow = plot3([-1.2 1.2],[0 0],[0 0]); 
h_x_arrow_a = plot3([1.15 1.2],[0.05 0],[0 0]); 
h_x_arrow_b = plot3([1.15 1.2],[-0.05 0],[0 0]); 

h_y_arrow = plot3([0 0],[-1.2 1.2],[0 0]); 
h_y_arrow_a = plot3([0.05 0],[1.15 1.2],[0 0]);
h_y_arrow_b = plot3([-0.05 0],[1.15 1.2],[0 0]);

h_z_arrow = plot3([0 0],[0 0],[0 1.2]); 
h_z_arrow_a = plot3([0 0],[0.05 0],[1.15 1.2]);
h_z_arrow_b = plot3([0 0],[-0.05 0],[1.15 1.2]);


h_ops_pt =  plot3([-0.25 -0.5],[0.25 0.5],[0 0.7]); 
h_c = patch(-0.5 + 0.01*cos(2*pi*s),0.5 + 0.01*sin(2*pi*s),0.7* ...
            ones(length(s),1),[1 1 1]);

xlabel('x');
ylabel('y');
zlabel('z');
