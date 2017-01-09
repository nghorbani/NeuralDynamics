% Nima Ghorbani
% Goals > HH Multicompartment active neurite

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% 
dt = 0.025;   %ms
te = 60;      %ms
ts = 260;     %ms
tend = 400;   %ms

%%%% Values from Book page 61
Cm = 1;    % uF
ENa = 50;  % mV
EK = -77;  % mV
EL = -54.4;% mV
gNa = 120; % mS
gK = 36;   % mS
gL = 0.3;  % mS
gax = 0.45; % mS

%%%% Values from Book page 61
alpha_m = @(v) (v==-40) + (~(v==-40))*eval_f(0.1*((v+40)/(1-exp(-(v+40)/10))));
beta_m = @(v) 4*exp(-(v+65)/18);

alpha_h = @(v) (0.07)*exp(-(v+65)/20);
beta_h = @(v) 1/(1+exp(-(v+35)/10));

alpha_n = @(v) (v==-55)*0.1 + (~(v==-55))*eval_f(0.01*((v+55)/(1-exp(-(v+55)/10))));
beta_n = @(v) 0.125*exp(-(v+65)/80);

m = @(v, m_prev) dt*[alpha_m(v)*(1-m_prev)-beta_m(v)*m_prev]+m_prev;
h = @(v, h_prev) dt*[alpha_h(v)*(1-h_prev)-beta_h(v)*h_prev]+h_prev;
n = @(v, n_prev) dt*[alpha_n(v)*(1-n_prev)-beta_n(v)*n_prev]+n_prev;
V = @(Ie,v_prev,v_prev_left,v_prev_right,m,h,n) (dt/Cm)*(Ie-gL*(v_prev-EL)-gNa*m^3*h*(v_prev-ENa)-gK*n^4*(v_prev-EK)-gax*(v_prev-v_prev_left)-gax*(v_prev-v_prev_right))+v_prev/Cm;

%%%%%%%%%%%%%%%

N = 100;   % Number of Compartments
je = 14;   % Ie injection compartment

[X Y] = meshgrid(dt:dt:tend,1:N);

I0s = [6 8 15 20]; % uA
%I0s = [8];

v =  zeros(N,tend/dt);
mt = zeros(N,tend/dt);
ht = zeros(N,tend/dt);
nt = zeros(N,tend/dt);

%setting correct initial values
%calculated from previous runs
v(1:N,1:te/dt) = -64.8525*ones(N,te/dt);
mt(1:N,1:te/dt) = 0.0529*ones(N,te/dt);
ht(1:N,1:te/dt) = 0.5954*ones(N,te/dt);
nt(1:N,1:te/dt) = 0.3177*ones(N,te/dt);
%%%


for figure_index=1:length(I0s)
    I0 = I0s(figure_index);
    Ie = zeros(N,tend/dt); % setup current matrix as zero at every time at ever compartment
    Ie(je,te/dt:(ts-te)/dt) = I0;
    for t=te:tend/dt % for every step of time  
        % calculate v,m,h,n at injection compartment
        v(je,t) = V(Ie(je,t),v(je,t-1),v(je-1,t-1),v(je+1,t-1),mt(je,t-1),ht(je,t-1),nt(je,t-1));
        mt(je,t) = m(v(je,t-1),mt(je,t-1));
        ht(je,t) = h(v(je,t-1),ht(je,t-1));
        nt(je,t) = n(v(je,t-1),nt(je,t-1));
        % then propagate
        % Forward propagation
        for j = je+1:N-1 % forward propagation 
            v(j,t) = V(Ie(j,t),v(j,t-1),v(j-1,t-1),v(j+1,t-1),mt(j,t-1),ht(j,t-1),nt(j,t-1));
            mt(j,t) = m(v(j,t-1),mt(j,t-1));
            ht(j,t) = h(v(j,t-1),ht(j,t-1));
            nt(j,t) = n(v(j,t-1),nt(j,t-1));
        end      
        % killed end boundary condition for last compartment
        % last compartment will always be grounded
        v(N,t) = -64.8525;% Em
        %Backward propagation
        for j = je-1:-1:2 % backward propagation 
            v(j,t) = V(Ie(j,t),v(j,t-1),v(j-1,t-1),v(j+1,t-1),mt(j,t-1),ht(j,t-1),nt(j,t-1));
            mt(j,t) = m(v(j,t-1),mt(j,t-1));
            ht(j,t) = h(v(j,t-1),ht(j,t-1));
            nt(j,t) = n(v(j,t-1),nt(j,t-1));
        end  
        % sealed end boundary for the first compartment
        % set equal to the right one
        v(1,t) = v(2,t);
    end
    figure(120+figure_index);
    surf(X,Y,v);shading interp;colormap hot;
    title(sprintf('Rectangular Impulse Input Current I_{0}= %d\\muA\nPotential Propagation In an Active Neurite',I0));

    xlabel('Time (ms)'); ylabel('Compartment (j)');zlabel('Voltage (mV)');
    ylim([1,N]);
end


