% Nima Ghorbani
% Goal > Analysis of linear dynamical DStems

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% Solving the DS
A = [-0.5,-0.5,0;-0.5,-0.5,0;0,0,2];

X0s = [[1,1,0];[100,0,0];[0,1,0];[0,0,1e-6]];
[Q,Lambda] = eig(A);
Xt = @(t,X0) Q*expm(Lambda*t)*Q^-1*X0';

ts = 0:0.1:10;
for i=1:4
    X0 = X0s(i,:);
    x = zeros(3,numel(ts));
    for j=1:numel(ts)
        x(:,j) =  Xt(ts(j),X0);
    end
    figure(100+i);grid on;
    subplot(3,2,1);plot(x(1,:),x(2,:));xlabel('X1');ylabel('X2');
    title(sprintf('1.1 X0 = %d %d %d',X0));
    subplot(3,2,3);plot(x(1,:),x(3,:));xlabel('X1');ylabel('X3');
    subplot(3,2,5);plot(x(2,:),x(3,:));xlabel('X2');ylabel('X3');
    subplot(3,2,[2,4,6]);plot3(x(1,:),x(2,:),x(3,:),'b');xlabel('X1');ylabel('X2');zlabel('X3');grid on;
    hold on;eq_point = plot3(x(1,end),x(2,end),x(3,end),'rO');
    legend(eq_point,sprintf('eq. point\n %1.2f %1.2f %1.2f',x(1,end),x(2,end),x(3,end)),'Location','northoutside');
end

%%
f = @(X) A*[X(1);X(2);X(3)];
t1 = linspace(-10,10,20);
t2 = linspace(-10,10,20);

[x,y] = meshgrid(t1,t2);

u = zeros(size(x));
v = zeros(size(x));

%computing derivative at each point
% x3 = 0
for i = 1:numel(x)
    Xp = f([x(i); y(i);0]);
    u(i) = Xp(1);
    v(i) = Xp(2);
end
figure(105);quiver(x,y,u,v,'r');
xlabel('x_1');ylabel('x_2');
axis tight equal;
title('Vector Field of Dynamics x_3=0');

% x2 = 0
for i = 1:numel(x)
    Xp = f([x(i); 0;y(i)]);
    u(i) = Xp(1);
    v(i) = Xp(3);
end
figure(106);quiver(x,y,u,v,'r');
xlabel('x_1');ylabel('x_3');
axis tight equal;
title('Vector Field of Dynamics x_2=0');


% x1 = 0
for i = 1:numel(x)
    Xp = f([0;x(i);y(i)]);
    u(i) = Xp(2);
    v(i) = Xp(3);
end
figure(107);quiver(x,y,u,v,'r');
xlabel('x_2');ylabel('x_3');
axis tight equal;
title('Vector Field of Dynamics x_1=0');

% whole phase plane

x1dom = linspace(-5,5,10);
x2dom = linspace(-5,5,10);
x3dom = linspace(-5,5,10);

[X1,X2,X3] = meshgrid(x1dom,x2dom,x3dom);

U = -.5*X1 -.5*X2; % dx1/dt
V = -.5*X1 -.5*X2; % dx2/dt
W = 2*X3;          % dx3/dt

figure(108);quiver3(X1,X2,X3,U,V,W,'r');axis tight equal;
xlabel('X1');ylabel('X2');zlabel('X3');
title('Vector Field of Whole Dynamics');
% autoArrangeFigures();


%% 

f = @(X) A*[X(1);X(2);X(3)];

e1=(1/sqrt(2))*[1 -1 0];
e2=[0 0 1];

[x,y] = meshgrid(linspace(-10,10,20),linspace(-10,10,20));

u = zeros(size(x));u_p = zeros(size(x));
v = zeros(size(x));v_p = zeros(size(x));

%computing derivative at each point
% also the projections
% x3=0
for i = 1:numel(x)
    Xp = f([x(i); y(i);0]);
    Xp_proj = f([x(i)*e1(1)/norm(e1); y(i)*e2(2)/norm(e2);0]);
    u(i) = Xp(1);
    v(i) = Xp(2);
    u_p(i) = Xp_proj(1)*e1(1)/norm(e1);
    v_p(i) = Xp_proj(2)*e2(2)/norm(e2);
end
figure(109);
subplot(1,2,1);quiver(x,y,u,v,'r');title('Vector Field of Dynamics x_3=0')
xlabel('x_1');ylabel('x_2');axis tight equal;

subplot(1,2,2);quiver(x,y,u_p,v_p,'r');title('Projected Vector Field of Dynamics x_3=0')
xlabel('Y_1');ylabel('Y_2');axis tight equal;

% x2=0
for i = 1:numel(x)
    Xp = f([x(i);0;y(i)]);
    Xp_proj = f([x(i)*e1(1)/norm(e1);0; y(i)*e2(3)/norm(e2)]);
    u(i) = Xp(1);v(i) = Xp(3);
    u_p(i) = Xp_proj(1)*e1(1)/norm(e1);
    v_p(i) = Xp_proj(3)*e2(3)/norm(e2);
end
figure(110);
subplot(1,2,1);quiver(x,y,u,v,'r');title('Vector Field of Dynamics x_2=0')
xlabel('x_1');ylabel('x_3');axis tight equal;

subplot(1,2,2);quiver(x,y,u_p,v_p,'r');title('Projected Vector Field of Dynamics x_2=0')
xlabel('Y_1');ylabel('Y_3');axis tight equal;

% x1 = 0
for i = 1:numel(x)
    Xp = f([0;x(i);y(i)]);
    Xp_proj = f([0; x(i)*e1(2)/norm(e1); y(i)*e2(3)/norm(e2)]);
    u(i) = Xp(2);v(i) = Xp(3);
    u_p(i) = Xp_proj(2)*e1(2)/norm(e1);
    v_p(i) = Xp_proj(3)*e2(3)/norm(e2);
end
figure(111);
subplot(1,2,1);quiver(x,y,u,v,'r');title('Vector Field of Dynamics x_1=0')
xlabel('x_2');ylabel('x_3');axis tight equal;

subplot(1,2,2);quiver(x,y,u_p,v_p,'r');title('Projected Vector Field of Dynamics x_1=0')
xlabel('Y_2');ylabel('Y_3');axis tight equal;

%%
A = [-0.5,-0.5,0;-0.5,-0.5,0;0,0,2];
S = [1 2 0]';

[Q,D] = eig(A);

f = @(X) A*[X(1);X(2);X(3)]+S;

f_proj = @(X) D*Q^-1*[X(1);X(2);X(3)]+Q^-1*S;


[x,y] = meshgrid(linspace(-10,10,20),linspace(-10,10,20));

X1 = zeros(size(x));Y1 = zeros(size(x));
X2 = zeros(size(x));Y2 = zeros(size(x));
X3 = zeros(size(x));Y3 = zeros(size(x));

%computing derivative at each point
% also the projections
% Y3=0
for i = 1:numel(x)
    Xp = f([x(i); y(i);0]);
    Yp = f_proj([x(i); y(i);0]);
    
    X1(i) = Xp(1);
    X2(i) = Xp(2);
    
    Y1(i) = Yp(1);
    Y2(i) = Yp(2);
    
end
figure(112);
subplot(1,2,1);quiver(x,y,X1,X2,'r');title('Vector Field of Dynamics X_3=0')
xlabel('x_1');ylabel('x_2');axis tight equal;

subplot(1,2,2);quiver(x,y,Y1,Y2,'r');title('Vector Field of Projected Dynamics Y_3=0')
xlabel('Y_1');ylabel('Y_2');axis tight equal;

%computing derivative at each point
% also the projections
% Y2=0
for i = 1:numel(x)
    Xp = f([x(i); 0;y(i)]);
    Yp = f_proj([x(i);0 ;y(i)]);

    X1(i) = Xp(1);
    X3(i) = Xp(3);
    
    Y1(i) = Yp(1);
    Y3(i) = Yp(3);

end
figure(113);
subplot(1,2,1);quiver(x,y,X1,X3,'r');title('Vector Field of Dynamics X_2=0')
xlabel('X_1');ylabel('X_3');axis tight equal;

subplot(1,2,2);quiver(x,y,Y1,Y3,'r');title('Vector Field of Projected Dynamics Y_2=0')
xlabel('Y_1');ylabel('Y_3');axis tight equal;

%computing derivative at each point
% also the projections
% Y1=0
for i = 1:numel(x)
    Xp = f([x(i); 0;y(i)]);
    Yp = f_proj([x(i);0 ;y(i)]);

    X2(i) = Xp(2);
    X3(i) = Xp(3);
    
    Y2(i) = Yp(2);
    Y3(i) = Yp(3);

end
figure(114);
subplot(1,2,1);quiver(x,y,X2,X3,'r');title('Vector Field of Dynamics X_1=0')
xlabel('X_1');ylabel('X_3');axis tight equal;

subplot(1,2,2);quiver(x,y,Y2,Y3,'r');title('Vector Field of Projected Dynamics Y_1=0')
xlabel('Y_1');ylabel('Y_3');axis tight equal;


%%
figure(115);hold on;
c = 2;s1 = 1;s2 = 1;
p_min = -5; p_max = 5; v_step = 0.5; p_step = 2.5;
DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);
xlabel('U1');ylabel('U2');
title('Vector field for a decision network')
plot(1/3,1/3,'bO',1,-1,'bO',-1,1,'bO');

%%
figure(116);hold on;
c = 2;s1 = 1.2;s2 = 1;
p_min = -5; p_max = 5; v_step = 0.2; p_step = 2.5;
DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);
plot(0.26,0.46,'bO',-.8,1,'bO',1.2,-1.4,'bO');

xlabel('U1');ylabel('U2');
title('Vector field for a DN s_1 = 1.2, s_2 = 1')

figure(126);hold on;
c = 2;s1 = 1;s2 = 1.2;
DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title('Vector field for a DN s_1 = 1, s_2 = 1.2')
%%

figure(117);hold on;
xlabel('U1');ylabel('U2');
title('Stationary Solutions for s_1 = 1.2, s_2 = 1')

c = 2;s1 = 1.2;s2 = 1;
p_min = -5; p_max = 5; v_step = 0.25; p_step = -1;

DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);
[ts,xs] = ode45(DS,[0 100],[0 0]);
plot(xs(:,1),xs(:,2),'b','LineWidth',1);p1=plot(xs(end,1),xs(end,2),'bO');
[ts,xs] = ode45(DS,[0 100],[1,-1]);
plot(xs(:,1),xs(:,2),'k','LineWidth',1);p2=plot(xs(end,1),xs(end,2),'kO');
[ts,xs] = ode45(DS,[0 100],[-1,1]);
plot(xs(:,1),xs(:,2),'g','LineWidth',1);p3=plot(xs(end,1),xs(end,2),'gO');
legend([p1 p2 p3],'U_0_,_1 = [0 0]','U_0_,_2 = [1 -1]','U_0_,_3 = [-1 1]');

figure(118);hold on;
xlabel('U1');ylabel('U2');
title('Stationary Solutions for s_1 = 1, s_2 = 1.2')

c = 2;s1 = 1;s2 = 1.2;

DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);
[ts,xs] = ode45(DS,[0 100],[0 0]);
plot(xs(:,1),xs(:,2),'b','LineWidth',1);p1=plot(xs(end,1),xs(end,2),'bO');
[ts,xs] = ode45(DS,[0 100],[1,-1]);
plot(xs(:,1),xs(:,2),'k','LineWidth',1);p2=plot(xs(end,1),xs(end,2),'kO');
[ts,xs] = ode45(DS,[0 100],[-1,1]);
plot(xs(:,1),xs(:,2),'g','LineWidth',1);p3=plot(xs(end,1),xs(end,2),'gO');
legend([p1 p2 p3],'U_0_,_1 = [0 0]','U_0_,_2 = [1 -1]','U_0_,_3 = [-1 1]');
%xs(end,:);

%%
figure(119);hold on;

c = -2;s1 = 1;s2 = 1;
p_min = -5; p_max = 5; v_step = .5; p_step = 2.5;
DS = @(t,x) [-x(1) - c*max(x(2),0)+s1;-x(2)-c*max(x(1),0)+s2];%t for trajectory plotting
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title('Vector field for a DN c = -2, s_1 = 1, s_2 = 1')
%% 
figure(120);hold on;

c = 2;s1 = 1;s2 = 1;
step_tresh = @(x) (x<=0)*0 + (x>0)*1;

p_min = -5; p_max = 5; v_step = .5; p_step = 2.5;
DS = @(t,x) [-x(1) - c*step_tresh(x(2))+s1;-x(2)-c*step_tresh(x(1))+s2];
vf_2d(DS,p_min,p_max,v_step,p_step);

xlabel('U1');ylabel('U2');
title('Vector field when using a step Treshold')
