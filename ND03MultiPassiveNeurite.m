% Nima Ghorbani
% Goal > Simulation of multi-compartment model

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% variables setup
dt = 1e-4; % delta t

Rm = 1.59e9;   % 1.59GOhm
Ra = 0.0318e9; % 0.0318GOhm
C = 62.8e-12;  % 62.8pF


%numerically claculating V at compartment j at time t with forward euler method
P = (C*Ra)^-1;
Q = ((dt^-1)-((C*Rm)^-1)-(2*((C*Ra)^-1)));

V = @(ie, v_previous, v_left_previous, v_rigth_previous) dt*(v_previous*Q + v_left_previous*P + v_rigth_previous*P + ie/C);

%%

%set timelapse
te = 20e-3; %seconds
ts = 1;     %seconds
tend = 1 ;  %seconds

% injected current
N = 50; %number of compartments
je = 20; % the compartment where we injected current

I0 = 10e-12; % -10pA
Ie = zeros(N,tend/dt); % setup current matrix as zero at every time at ever compartment
Ie(je,te/dt:(ts-te)/dt) = I0;

%compartments = [je ]% [je forward_propagate backward_propagate]

vt = zeros(N,tend/dt); % v timecourse at compartment j at time t
for t= 2:tend/dt % for every time step
    %calculate v(j,t) at injection compartment
    vt(je,t) = V(Ie(je,t), vt(je,t-1), vt(je-1,t-1), vt(je+1,t-1));
    %forward propagation
    for j = je:N-1 %last compartment will always be grounded
        vt(j,t) = V(Ie(j,t), vt(j,t-1), vt(j-1,t-1), vt(j+1,t-1));
    end
    %backward propagation
    for j = je:-1:2 %special treatment for first compartment
        vt(j,t) = V(Ie(j,t), vt(j,t-1), vt(j-1,t-1), vt(j+1,t-1));
    end
    vt(1,t) = vt(2,t);
end
figure(100);
[X Y] = meshgrid(dt:dt:tend,1:N);
surf(X,Y,vt);shading interp;colormap hot;
title(sprintf('Step Impulse Input Current\nPotential Propagation In a Passive Neurite'));
xlabel('Time (s)'); ylabel('Compartment (j)');zlabel('Voltage (V)');
ylim([1,N]);


%%

%set timelapse
te = 20e-3; %seconds
ts = 400e-3;     %seconds
tend = 1;  %seconds

% injected current
N = 50; %number of compartments
je = 14; % the compartment where we injected current

I0 = 100e-12; % -10pA
Ie = zeros(N,tend/dt); % setup current matrix as zero at every time at ever compartment
Ie(je,te/dt:(ts-te)/dt) = I0;

%compartments = [je ]% [je forward_propagate backward_propagate]

vt = zeros(N,tend/dt); % v timecourse at compartment j at time t
for t= 2:tend/dt % for every time step
    %calculate v(j,t) at injection compartment
    vt(je,t) = V(Ie(je,t), vt(je,t-1), vt(je-1,t-1), vt(je+1,t-1));
    %forward propagation
    for j = je:N-1 %last compartment will always be grounded
        vt(j,t) = V(Ie(j,t), vt(j,t-1), vt(j-1,t-1), vt(j+1,t-1));
    end
    %backward propagation
    for j = je:-1:2 %special treatment for first compartment
        vt(j,t) = V(Ie(j,t), vt(j,t-1), vt(j-1,t-1), vt(j+1,t-1));
    end
    vt(1,t) = vt(2,t);
end
figure(101);
surf(X,Y,vt);shading interp;colormap hot;
title(sprintf('Rectangular Impulse Input Current\nPotential Propagation In a Passive Neurite'));
xlabel('Time (s)'); ylabel('Compartment (j)');zlabel('Voltage (V)');
ylim([1,N]);