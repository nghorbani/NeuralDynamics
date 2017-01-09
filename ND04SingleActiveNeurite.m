% Nima Ghorbani
% Goal > Sinle Active Neuron

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');


%% Values from Principles of Computational Mod - David Sterratt
clear;

dt = 0.025;   %ms
te = 50;      %ms
ts = 300;     %ms
tend = 1000;  %ms

%%%%%% page 61
Cm = 1;    % uF
ENa = 50;  % mV
EK = -77;  % mV
EL = -54.4;% mV
gNa = 120; % mS
gK = 36;   % mS
gL = 0.3;  % mS


%%%%%% page 61
alpha_m = @(v) (v==-40) + (~(v==-40))*eval_f(0.1*((v+40)/(1-exp(-(v+40)/10))));
beta_m = @(v) 4*exp(-(v+65)/18);

alpha_h = @(v) (0.07)*exp(-(v+65)/20);
beta_h = @(v) 1/(1+exp(-(v+35)/10));

alpha_n = @(v) (v==-55)*0.1 + (~(v==-55))*eval_f(0.01*((v+55)/(1-exp(-(v+55)/10))));
beta_n = @(v) 0.125*exp(-(v+65)/80);


m = @(v, m_prev) dt*[alpha_m(v)*(1-m_prev)-beta_m(v)*m_prev]+m_prev;
h = @(v, h_prev) dt*[alpha_h(v)*(1-h_prev)-beta_h(v)*h_prev]+h_prev;
n = @(v, n_prev) dt*[alpha_n(v)*(1-n_prev)-beta_n(v)*n_prev]+n_prev;
V = @(Ie,v_prev,m,h,n) (dt/Cm)*(Ie-gL*(v_prev-EL)-gNa*m^3*h*(v_prev-ENa)-gK*n^4*(v_prev-EK))+v_prev/Cm;

%%%debug
figure(102);hold on;

v_temp = -90:0.0001:100; %mV
m_temp = zeros(1,length(v_temp));
h_temp = zeros(1,length(v_temp));
n_temp = zeros(1,length(v_temp));
for t=2:length(v_temp)
    m_temp(t) = m(v_temp(t-1),m_temp(t-1));
    h_temp(t) = h(v_temp(t-1),h_temp(t-1));
    n_temp(t) = n(v_temp(t-1),n_temp(t-1));
end

plot(v_temp,m_temp,'-');plot(v_temp,h_temp,'.');plot(v_temp,n_temp,'--');
legend('m','h','n');
ylabel('Probability'); xlabel('V_{mV}');
xlim([min(v_temp),max(v_temp)]);
title('Steady State Ion Channel Gating Probability')
%%%%%%%%%%%%%%%

I0s = [6]; % uA

v =  zeros(length(I0s),tend/dt);
mt = zeros(length(I0s),tend/dt);
ht = zeros(length(I0s),tend/dt);
nt = zeros(length(I0s),tend/dt);

%setting correct initial values
%calculated from previous runs
v(1:length(I0s),1:te/dt) = -64.8525*ones(length(I0s),te/dt);
mt(1:length(I0s),1:te/dt) = 0.0529*ones(length(I0s),te/dt);
ht(1:length(I0s),1:te/dt) = 0.5954*ones(length(I0s),te/dt);
nt(1:length(I0s),1:te/dt) = 0.3177*ones(length(I0s),te/dt);
%%%

figure(103);

for k=1:length(I0s)
    I0 = I0s(k);    
    Ie = zeros(1,tend/dt); % setup current matrix as zero at every time at ever compartment
    Ie(1,te/dt:(ts-te)/dt) = I0;
    
    for t=2:tend/dt % v(0) = 0        
        v(k,t) = V(Ie(t),v(k,t-1),mt(k,t-1),ht(k,t-1),nt(k,t-1));
        mt(k,t) = m(v(k,t-1),mt(k,t-1));
        ht(k,t) = h(v(k,t-1),ht(k,t-1));
        nt(k,t) = n(v(k,t-1),nt(k,t-1));
    end      
    subplot(length(I0s),1,k);plot(dt:dt:tend,v(k,:));
    ylabel('V_{mV}');
    legend(sprintf('I_{0} = %d \\muA',I0));
end
xlabel('t_{ms}'); 
subplot(length(I0s),1,1);title('HH Model');

%% 

dt = 0.025;  %ms
te = 50;     %ms
ts = 800;    %ms
tend = 1000; %ms

I0s = [0:0.5:20]; % uA

v = zeros(length(I0s),tend/dt);
mt = zeros(length(I0s),tend/dt);
ht = zeros(length(I0s),tend/dt);
nt = zeros(length(I0s),tend/dt);
rates = zeros(length(I0s),1);

T = (ts-te)*10^-3; % ms
for k=1:length(I0s)
    I0 = I0s(k);
    Ie = [zeros(1,te/dt) I0*ones(1,uint32((ts-te)/dt)) zeros(1,uint32((tend-ts)/dt))];
    for t=2:tend/dt % v(0) = 0        
        v(k,t) = V(Ie(t),v(k,t-1),mt(k,t-1),ht(k,t-1),nt(k,t-1));
        mt(k,t) = m(v(k,t-1),mt(k,t-1));
        ht(k,t) = h(v(k,t-1),ht(k,t-1));
        nt(k,t) = n(v(k,t-1),nt(k,t-1));
    end     
    rates(k) = find_rate(v(k,:),20,T) ;%threshold 20mV
end
figure(104);
plot(I0s,rates);
xlabel('I_{0\muA}'); ylabel('Rate (Spikes/ms)');
title('Firing rate for various I_0')
