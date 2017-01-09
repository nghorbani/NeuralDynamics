% Nima Ghorbani
% Goal > Simulation of two-compartment model

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% variables setup
dt = 1e-3; % delta t
em = 0; % membrane equilibrium potential, v

Rm = 265e6;  % 265MOhm
C = 75e-12; % 75 pF
Ras = [7e6 265e6 30e9];    % 7 MOhm
 
te = 0.4; %seconds
ts = 0.44; %seconds
tend = 1 ; %seconds

%determining timecourse for V(t)


%numerically claculating v1 and v2
v1 = @(ie, P, K, Ra, v1_previous, v2_previous) (1/P)*((v1_previous/dt)+(((1/K)*v2_previous)/(C*Ra*dt))+(ie/C));%vt numerical form with backward euler method
v2 = @(ie, K, Ra, v1, v2_previous) (1/K)*((v1/(C*Ra) + v2_previous/dt));%vt numerical form with backward euler method

%% 3.2

% injected current
I0 = -100e-12; % -100pA
Ie = [zeros(1,te/dt) I0*ones(1,uint8((ts-te)/dt)) zeros(1,(tend-ts)/dt)];

figure(311); % figure for v1

str_legend =  {}; % array for legends

for figure_i=1:length(Ras)
    Ra = Ras(figure_i);
    K = (1/dt)+(1/(C*Rm))+(1/(C*Ra));
    P = K - (1/(K*C*C*Ra*Ra));
    %setup t and Ie
    v1_t = zeros(1,tend/dt);% v1
    v2_t = zeros(1,tend/dt);% v2
    for t=2:tend/dt %inital v is zero, so counting from 2
        v1_t(t) = v1(Ie(t),P, K, Ra, v1_t(t-1), v2_t(t-1));
        v2_t(t) = v2(Ie(t), K, Ra, v1_t(t), v2_t(t-1));
    end

    subplot(1,3,figure_i),hold on;
    plot(dt:dt:tend,v1_t);plot(dt:dt:tend,v2_t);
    legend('v1','v2');
    xlabel('time (s)');ylabel('amplitude (v)');
    title(sprintf('Ra = %d M\\Omega',Ra*(10^-6)));

end

%% 3.3
Ra = 300e6;  % 300MOhm

dt = 1e-4; % changing delta for high frequencies

K = (1/dt)+(1/(C*Rm))+(1/(C*Ra));
P = K - (1/(K*C*C*Ra*Ra));

tend = 5;

freqs = [1 2 5 10 20 50 100 200 500 1e3 2e3 5e3]; % frequencies in Herz

I0 = 100e-12;
Ie = I0*sin(2*pi*freqs'*(0:dt:tend));

for freq_i=1:length(freqs)

    %%%%%%%%%%%
    v1_t = zeros(1,tend/dt);% v1
    v2_t = zeros(1,tend/dt);% v2
    for t=2:tend/dt %inital v is zero, so counting from 2
        v1_t(t) = v1(Ie(freq_i,t), P, K, Ra, v1_t(t-1), v2_t(t-1));
        v2_t(t) = v2(Ie(freq_i,t), K, Ra, v1_t(t), v2_t(t-1));
    end

    freq_responses(1,freq_i) = max(v1_t(1,1000:end));
    freq_responses(2,freq_i) = max(v2_t(1,1000:end));
    
end
figure(331); hold on;
loglog(freqs, freq_responses(1,:),'r');%plot(log(freqs), log(freq_responses(1,:)),'*r');


loglog(freqs, freq_responses(2,:),'b');%plot(log(freqs), log(freq_responses(2,:)),'*b');


ylabel('Log Amplitude');
xlabel('Log Frequency');
title('3.3 - Bode Diagram for 3 M\Omega ')
legend('V1', 'V2');

