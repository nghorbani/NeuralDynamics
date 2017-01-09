% Nima Ghorbani
% Goal > Simulation of single compartment model

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% variables setup
dts = [1e-4 1e-3 10e-3]; % timesteps, s
em = 0; % membrane equilibrium potential, v
l = 100e-6; % length, m 
d = 2e-6; % diameter, m
r_m = 1; % membrane specific resistance, ohm*m^2
r_a = 1; % axial specific resistance, ohm*m
c_m = 1e-2; % membrane specific capaitance, F*m^-2

R_m = r_m/(pi*d);%(4*pi*(d/2)^2); % ohm, Resistance = specific_r / surface area;
C_m = c_m *(pi*d); %(4*pi*(d/2)^2); % farad, capacitance = specific_c * surface area;
 
start_time = 0;
end_time = 5;

%% 
%determining timecourse for V(t)
%plotting for different dt



%claculating vt
v = @(dt, R, C, ie, em, vt_previous) (R/(R*C + dt))*(ie*dt + vt_previous*C);%vt numerical form with backward euler method
% injected current
I_inject = @(time_course) (-50e-12)*ones(1,length(time_course));% 50 pico ampere for t>0


figure(100),hold on; % create figure
str_legend =  {}; % array for legends

for dt=dts
    %setup t and Ie
    time_course = start_time:dt:end_time;%%setup time
    Ie = I_inject(time_course);% 50 pico ampere for t>0
    %%%%%%%%%%%%
    vt = zeros(1,length(time_course));% vt
    for i=2:length(time_course)%inital v is zero, so counting from 2
        vt(i) = v(dt, R_m, C_m, Ie(i), em, vt(i-1));
    end

    plot(time_course,vt);
    str_legend = [str_legend, sprintf('dt = %2.2fms',dt*1000.)];
    
end
legend(str_legend);
xlabel('time');
ylabel('voltage');
%title('Excercises 3.1 & 3.2')
%axis([-dt 0.3 0 -4]); % set axis for better visualization

%% 
dt = 0.1e-3; % dt = 0.1 ms

figure(200),hold on; % create figure

%%%%% new specific membrane capacitance

c_m = 1e-1; % new membrane specific capaitance
C_m = c_m *(pi*d); %(4*pi*(d/2)^2); % farad, capacitance = specific_c * surface area;


%setup t and Ie
time_course = start_time:dt:end_time;%%setup time
Ie = I_inject(time_course);% -50 pico ampere for t>0

vt = zeros(1,length(time_course));% vt
for i=2:length(time_course)%inital v is zero, so counting from 2
    vt(i) = v(dt, R_m, C_m, Ie(i), em, vt(i-1));
end

plot(time_course,vt);

%%%%%% new specific membrane resistance


r_m = 10; % membrane specific resistance, ohm*m^2
R_m = r_m/(pi*d);%(4*pi*(d/2)^2); % ohm, Resistance = specific_r / surface area;

%setup t and Ie
time_course = start_time:dt:end_time;%%setup time
Ie = I_inject(time_course);% -50 pico ampere for t>0
%%%%%%%%%%%%
vt = zeros(1,length(time_course));% vt
for i=2:length(time_course)%inital v is zero, so counting from 2
    vt(i) = v(dt, R_m, C_m, Ie(i), em, vt(i-1));
end

plot(time_course,vt);

%%%%%% plot values
%title('Excercise 3.3')
legend('c_{m} = 10^-^1 Fm^-^2','r_{m} = 10 \Omega*m^2');
xlabel('time');
ylabel('voltage');
%axis([-dt 0.3 0 4]); % set axis for better visualization

%% sinusoidal input and bode diagram
figure(300),hold on; % create figure
str_legend =  {}; % array for legends

start_time = 0;
end_time = 5;

r_m = 1; % membrane specific resistance, ohm*m^2
c_m = 1e-2; % membrane specific capaitance, F*m^-2

R_m = r_m/(pi*d);%(4*pi*(d/2)^2); % ohm, Resistance = specific_r / surface area;
C_m = c_m *(pi*d); %(4*pi*(d/2)^2); % farad, capacitance = specific_c * surface area;

freqs = [.5, 1, 2, 8, 100, 1000 ]; % frequencies in Herz
dt = 0.1e-3; % dt = 0.1 ms
freq_responses = [];

for f=freqs
    %setup t and Ie
    time_course = start_time:dt:end_time;%%setup time
    Ie = (100e-12)*sin(2*pi*f*time_course);% -50 pico ampere for t>0
    %%%%%%%%%%%%
    vt = zeros(1,length(time_course));% vt
    for i=2:length(time_course)%inital v is zero, so counting from 2
        vt(i) = v(dt, R_m, C_m, Ie(i), em, vt(i-1));
    end

    plot(time_course,vt);

    str_legend = [str_legend, sprintf('f = %2.1fHz',f)];
    freq_responses = [freq_responses, max(round(abs(vt),2,'significant'))];
    
    
end
legend(str_legend);
xlabel('time');
ylabel('voltage');
%title('Excercise 3.4')

figure(400); % create figure
loglog(freqs, freq_responses);
ylabel('Amplitude of Response');
xlabel('Frequency');
%title('Excercise 3.4')

