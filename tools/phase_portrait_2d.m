% goal is to plot a 3d phase portrain for a 3*3 coeficient DS matrix
clear all; close all; clc;
%set by the desired domain
x1dom = linspace(-100,100,20);
x2dom = linspace(-100,100,20);


[X1,X2] = meshgrid(x1dom,x2dom); % generate mesh of domain

c = 2;s1 = 1;s2 = 1;

U = -X1 - c*max(X2,0)+s1; % dx1/dt
V = -X2 - c*max(X1,0)+s2; % dx2/dt

h = quiver(X1,X2,U,V);
set(h,'MaxHeadSize',1e5,'AutoScaleFactor',3);
axis tight equal;
xlabel('X1');
ylabel('X2');
