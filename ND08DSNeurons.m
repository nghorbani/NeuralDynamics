% Nima Ghorbani
% Goal > phase portraits for a nonlinear network with two divisive inhibitor neurons

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% Nonlinear network with two divisive inhibition neurons
figure(100);hold on;
tav = 0.025; s1 = 0; s2 = 0;
p_min = -2; p_max = 2; v_step = 0.2; p_step = 100;
DS = @(t,x) [(1/tav)*(-x(1)+(s1/(1+x(2))));(1/tav)*(-x(2)+(s1/(1+x(1))))];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title(sprintf('VF/PP for network of divisive inhibition neurons\ns = 0'))
plot(0,0,'bO');

figure(101);hold on;
tav = 0.025; s1 = 3/4; s2 = 3/4;
DS = @(t,x) [(1/tav)*(-x(1)+(s1/(1+x(2))));(1/tav)*(-x(2)+(s1/(1+x(1))))];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title('s = 3/4')
plot(0.5,0.5,'bO');

%% Simple autoassociative memory
figure(102);hold on;
tav = 0.025; 
p_min = -2; p_max = 2; v_step = 0.1; p_step = 0.5;
DS = @(t,u) [(1/tav)*(-u(1)+max(-.1*u(1)+u(2),0));(1/tav)*(-u(2)+max(-.1*u(1)-.1*u(2),0))];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U2');ylabel('U3');
title('U1=0')

figure(103);hold on;
DS = @(t,u) [(1/tav)*(-u(1)+max(u(1)-.1*u(2),0));(1/tav)*(-u(2)+max(-.1*u(1)-.1*u(2),0))];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U3');
title('U2=0')

figure(104);hold on;
DS = @(t,u) [(1/tav)*(-u(1)+max(u(1)-0.1*u(2),0));(1/tav)*(-u(2)+max(-.1*u(1)+u(2),0))];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title('U3=0')

%% 3D quiver plot
precision = 10;% steps
dom_min = -5; dom_max = 5;% setting computations domain size

x1dom = linspace(dom_min,dom_max,precision);
x2dom = linspace(dom_min,dom_max,precision);
x3dom = linspace(dom_min,dom_max,precision);

[X1,X2,X3] = meshgrid(x1dom,x2dom,x3dom); % generate mesh of domain

dX1 = -X1 + max(- X2/10 - X3/10,0); % dx1/dt
dX2 = -X2 + max(- X1/10 - X3/10,0); % dx2/dt
dX3 = -X3 + max(- X1/10 - X2/10 - (11*X3)/10,0);  %dx3/dt

figure(105);
h = quiver3(X1,X2,X3,dX1,dX2,dX3);
set(h,'MaxHeadSize',1e5,'AutoScaleFactor',2);
xlabel('X1');ylabel('X2');zlabel('X3');
axis tight equal;

% 
figure(106);hold on;
v2=-10:0.1:10;
for k = 0:10
    plot(v2,k+v2.^2/2-log(-1*v2));
end
xlabel('V1');ylabel('V2');
title('Countour Lines describing trajectories of the dynamics');