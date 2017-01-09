% Nima Ghorbani
% Goals > 1-HH model extended with with A-Type current 
%         2-HH Extended Gatings
%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% HH extended model
dt = 0.025;   %ms
te = 60;      %ms
ts = 460;     %ms
tend = 1000;  %ms

Cm = 1;    % uF
ENa = 50;  % mV
EK = -77;  % mV
EA = -80;  % mV
EL = -22;  % mV
V0 = -73;  % mV
gNa = 120; % mS
gK = 20;   % mS
gA = 47.7; % mS
gL = 0.3;  % mS

alpha_m = @(v) (v==-34.7)*3.8 + (~(v==-34.7))*eval_f(3.8*.1*((v+34.7)/(1-exp(-(v+34.7)/10))));
beta_m = @(v) 3.8*4*exp(-(v+59.7)/18);

alpha_h = @(v) 3.8*0.07*exp(-(v+53)/20);
beta_h = @(v) 3.8/(1+exp(-(v+23)/10));

alpha_n = @(v) (v==-50.7)*0.1*1.9 + (~(v==-50.7))*eval_f(1.9*0.01*((v+50.7)/(1-exp(-(v+50.7)/10))));
beta_n = @(v) 1.9*0.125*exp(-(v+60.7)/80);

a_inf = @(v) ((0.0761*exp((v+99.22)/(31.84)))/(1+exp((v+6.17)/28.93)))^(1/3);
T_a = @(v) 0.3632+1.158/(1+exp((v+60.96)/20.12));
b_inf = @(v) 1/(1+exp((v+58.3)/14.54))^4;
T_b = @(v) 1.24 + 2.678/(1+exp((v-55)/16.027));


m = @(v, m_prev) dt*[alpha_m(v)*(1-m_prev)-beta_m(v)*m_prev]+m_prev;
h = @(v, h_prev) dt*[alpha_h(v)*(1-h_prev)-beta_h(v)*h_prev]+h_prev;
n = @(v, n_prev) dt*[alpha_n(v)*(1-n_prev)-beta_n(v)*n_prev]+n_prev;
a = @(v, a_prev) dt*[(a_inf(v) - a_prev)/T_a(v)] + a_prev;
b = @(v, b_prev) dt*[(b_inf(v) - b_prev)/T_b(v)] + b_prev;
V = @(Ie,v_prev,a,b,m,h,n) (dt/Cm)*(Ie-gL*(v_prev-EL)-gNa*m^3*h*(v_prev-ENa)-gK*n^4*(v_prev-EK)-gA*a^3*b*(v_prev-EA))+v_prev/Cm;


I0s = [0:0.5:20]; % uA

v = zeros(length(I0s),tend/dt);
mt = zeros(length(I0s),tend/dt);
ht = zeros(length(I0s),tend/dt);
nt = zeros(length(I0s),tend/dt);
at = zeros(length(I0s),tend/dt);
bt = zeros(length(I0s),tend/dt);

N = length(I0s);
%setting correct initial values
%calculated from previous runs
v(1:N,1:te/dt) = -72.9689*ones(N,te/dt);
mt(1:N,1:te/dt) = 0.0101*ones(N,te/dt);
ht(1:N,1:te/dt) = 0.9659*ones(N,te/dt);
nt(1:N,1:te/dt) = 0.1560*ones(N,te/dt);
at(1:N,1:te/dt) = 0.5405*ones(N,te/dt);
bt(1:N,1:te/dt) = 0.2882*ones(N,te/dt);
%%%

rates = zeros(length(I0s),1);

T = (ts-te)*10^-3; % ms
for k=1:length(I0s)
    I0 = I0s(k);
    Ie = [zeros(1,te/dt) I0*ones(1,uint32((ts-te)/dt)) zeros(1,uint32((tend-ts)/dt))];
    for t=te:tend/dt % v(0) = 0  
        v(k,t) = V(Ie(t),v(k,t-1),at(k,t-1),bt(k,t-1),mt(k,t-1),ht(k,t-1),nt(k,t-1));
        mt(k,t) = m(v(k,t-1),mt(k,t-1));
        ht(k,t) = h(v(k,t-1),ht(k,t-1));
        nt(k,t) = n(v(k,t-1),nt(k,t-1));
        at(k,t) = a(v(k,t-1),at(k,t-1));
        bt(k,t) = b(v(k,t-1),bt(k,t-1));
    end     
    rates(k) = find_rate(v(k,:),10,T); % threshold 10mV
end
figure(102);
plot(I0s,rates);
xlabel('I_{0\muA}'); ylabel('Rate (Spikes/ms)');
title('Firing rate for various I0')

%% HH Extended Gatings
clear all;

dt = 0.025;   %ms
te = 60;      %ms
ts = 460;     %ms
tend = 1000;  %ms

Cm = 1;    % uF
ENa = 50;  % mV
EK = -77;  % mV
EA = -80;  % mV
EL = -22;  % mV
V0 = -73;  % mV
gNa = 120; % mS
gK = 20;   % mS
gA = 47.7; % mS
gL = 0.3;  % mS

alpha_m = @(v) (v==-34.7)*3.8 + (~(v==-34.7))*eval_f(3.8*.1*((v+34.7)/(1-exp(-(v+34.7)/10))));
beta_m = @(v) 3.8*4*exp(-(v+59.7)/18);

alpha_h = @(v) 3.8*0.07*exp(-(v+53)/20);
beta_h = @(v) 3.8/(1+exp(-(v+23)/10));

alpha_n = @(v) (v==-50.7)*0.1*1.9 + (~(v==-50.7))*eval_f(1.9*0.01*((v+50.7)/(1-exp(-(v+50.7)/10))));
beta_n = @(v) 1.9*0.125*exp(-(v+60.7)/80);

a_inf = @(v) ((0.0761*exp((v+99.22)/(31.84)))/(1+exp((v+6.17)/28.93)))^(1/3);
T_a = @(v) 0.3632+1.158/(1+exp((v+60.96)/20.12));
b_inf = @(v) 1/(1+exp((v+58.3)/14.54))^4;
T_b = @(v) 1.24 + 2.678/(1+exp((v-55)/16.027));


m = @(v, m_prev) dt*[alpha_m(v)*(1-m_prev)-beta_m(v)*m_prev]+m_prev;
h = @(v, h_prev) dt*[alpha_h(v)*(1-h_prev)-beta_h(v)*h_prev]+h_prev;
n = @(v, n_prev) dt*[alpha_n(v)*(1-n_prev)-beta_n(v)*n_prev]+n_prev;
a = @(v, a_prev) dt*[(a_inf(v) - a_prev)/T_a(v)] + a_prev;
b = @(v, b_prev) dt*[(b_inf(v) - b_prev)/T_b(v)] + b_prev;
V = @(Ie,v_prev,a,b,m,h,n) (dt/Cm)*(Ie-gL*(v_prev-EL)-gNa*m^3*h*(v_prev-ENa)-gK*n^4*(v_prev-EK)-gA*a^3*b*(v_prev-EA))+v_prev/Cm;


v = zeros(1,tend/dt);
mt = zeros(1,tend/dt);
ht = zeros(1,tend/dt);
nt = zeros(1,tend/dt);
at = zeros(1,tend/dt);
bt = zeros(1,tend/dt);

%setting correct initial values
%calculated from previous runs
v(1,1:te/dt) = -72.9689*ones(1,te/dt);
mt(1,1:te/dt) = 0.0101*ones(1,te/dt);
ht(1,1:te/dt) = 0.9659*ones(1,te/dt);
nt(1,1:te/dt) = 0.1560*ones(1,te/dt);
at(1,1:te/dt) = 0.5405*ones(1,te/dt);
bt(1,1:te/dt) = 0.2882*ones(1,te/dt);
%%%


% evaluation of the model
I0 = 9;%uA    
Ie = [zeros(1,te/dt) I0*ones(1,uint32((ts-te)/dt)) zeros(1,uint32((tend-ts)/dt))];

for t=te:tend/dt % v(0) = 0    
    % HH Extended model calculations
    v(1,t) = V(Ie(t),v(1,t-1),at(1,t-1),bt(1,t-1),mt(1,t-1),ht(1,t-1),nt(1,t-1));
    mt(1,t) = m(v(1,t-1),mt(1,t-1));
    ht(1,t) = h(v(1,t-1),ht(1,t-1));
    nt(1,t) = n(v(1,t-1),nt(1,t-1));
    at(1,t) = a(v(1,t-1),at(1,t-1));
    bt(1,t) = b(v(1,t-1),bt(1,t-1));
end  
% plot
figure(104);
subplot(2,1,1);plot(dt:dt:tend,v(1,:));
ylabel('V_{mV}');xlabel('t_{ms}'); 
legend(sprintf('I_{0} = %d \\muA',I0));
title('HH Extended Model with "A" Current')
%axis([90,160,-80,60]);
subplot(2,1,2);hold on;
plot(dt:dt:tend,mt(1,:),'r.');plot(dt:dt:tend,ht(1,:),'b--');plot(dt:dt:tend,nt(1,:),'k');
plot(dt:dt:tend,at(1,:),'g');plot(dt:dt:tend,bt(1,:),'y');

ylabel('Gating Probability');xlabel('t_{ms}'); 
legend('m- Na activation','h- Na inactivation ','n- K activation','a','b');
title('Variations in Gating Variables for HH Extended model')
%axis([90,160,0,1]);
