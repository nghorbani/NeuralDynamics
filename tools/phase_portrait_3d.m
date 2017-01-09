% goal is to plot a 3d phase portrain for a 3*3 coeficient DS matrix
% 3D quiver plot for dynamical systems

precision = 5;% steps
dom_min = -5; dom_max = 5;% setting computations domain size

x1dom = linspace(dom_min,dom_max,precision);
x2dom = linspace(dom_min,dom_max,precision);
x3dom = linspace(dom_min,dom_max,precision);

[X1,X2,X3] = meshgrid(x1dom,x2dom,x3dom); % generate mesh of domain

dX1 = -X1 + max(- X2/10 - X3/10,0); % dx1/dt
dX2 = -X2 + max(- X1/10 - X3/10,0); % dx2/dt
dX3 = -X3 + max(- X1/10 - X2/10 - (11*X3)/10,0);  %dx3/dt

h = quiver3(X1,X2,X3,dX1,dX2,dX3);
set(h,'MaxHeadSize',1e5,'AutoScaleFactor',2);
axis tight equal;
xlabel('X1');ylabel('X2');zlabel('X3');