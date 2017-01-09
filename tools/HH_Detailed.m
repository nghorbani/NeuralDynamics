% Goals: 
% 1 - Variability of different Gating Probabilities in time course of AP in HH
% model
% 2 - See the effect of Depolarizing Pre Pulse 

clear; close all;clc;
set(0,'DefaultFigureWindowStyle','normal');

% Stimulation Wave Shape: 0_______te--tdpp=======ts________tend
dt = 0.001;       %ms
te = 100;         %ms step input start
tdpp = te + 0;     %ms Depolarizing Pre Pulse time
ts = tdpp + 15;    %ms step input reset
tend = 200;       %ms

% Stimulation amplitude
I_amp = 9;     %uA  
Idpp = 0;   %Depolarizing Pre Pulse uA


%%%%%% Values from Book page 61
Cm = 1;    % uF
ENa = 50;  % mV
EK = -77;  % mV
EL = -54.4;% mV
gNa = 120; % mS
gK = 36;   % mS
gL = 0.3;  % mS

%%%%%% Values from Book page 61
alpha_m = @(v) (v==-40) + (~(v==-40))*eval_f(0.1*((v+40)/(1-exp(-(v+40)/10))));
beta_m = @(v) 4*exp(-(v+65)/18);

alpha_h = @(v) (0.07)*exp(-(v+65)/20);
beta_h = @(v) 1/(1+exp(-(v+35)/10));

alpha_n = @(v) (v==-55)*0.1 + (~(v==-55))*eval_f(0.01*((v+55)/(1-exp(-(v+55)/10))));
beta_n = @(v) 0.125*exp(-(v+65)/80);

m = @(v, m_prev) dt*[alpha_m(v)*(1-m_prev)-beta_m(v)*m_prev]+m_prev;
h = @(v, h_prev) dt*[alpha_h(v)*(1-h_prev)-beta_h(v)*h_prev]+h_prev;
n = @(v, n_prev) dt*[alpha_n(v)*(1-n_prev)-beta_n(v)*n_prev]+n_prev;
gna = @(m,h) gNa*m^3*h;
gk = @(n) gK*n^4;
V = @(Ie,v_prev,m,h,n) (dt/Cm)*(Ie-gL*(v_prev-EL)-gna(m,h)*(v_prev-ENa)-gk(n)*(v_prev-EK))+v_prev/Cm;

Ina = @(v,m,h) gna(m,h)*(v-ENa);
Ik  = @(v,n) gk(n)*(v-EK);
%%
v  = zeros(1,tend/dt);
mt = zeros(1,tend/dt); 
ht = zeros(1,tend/dt);
nt = zeros(1,tend/dt);
gNat = zeros(1,tend/dt); 
gKt = zeros(1,tend/dt);
Ilt = zeros(1,tend/dt); 
Inat = zeros(1,tend/dt); 
Ikt = zeros(1,tend/dt);

%setting correct initial values
%calculated from previous runs
v(1,1:te/dt) = -64.8525*ones(1,te/dt);
mt(1,1:te/dt) = 0.0529*ones(1,te/dt);
ht(1,1:te/dt) = 0.5954*ones(1,te/dt);
nt(1,1:te/dt) = 0.3177*ones(1,te/dt);
% % %%%
% % %%
% % 
% % % evaluation of the model
% % Ie = [zeros(1,te/dt) Idpp*ones(1,uint32((tdpp-te)/dt)) I_amp*ones(1,uint32((ts-tdpp)/dt)) zeros(1,uint32((tend-ts)/dt))];
% % 
% % for t=te:tend/dt % v(0) = 0    
% %     % HH model calculations
% %     v(1,t) = V(Ie(t),v(1,t-1),mt(1,t-1),ht(1,t-1),nt(1,t-1));
% %     
% %     % extra plots
% %     Inat(1,t) = Ina(v(1,t-1),mt(1,t-1),ht(1,t-1));
% %     Ikt(1,t)  = Ik(v(1,t-1),nt(1,t-1));
% %     gNat(1,t) = gna(mt(1,t-1),ht(1,t-1)); 
% %     gKt(1,t)  = gk(nt(1,t-1));
% %     % extra plots done
% % 
% %     mt(1,t) = m(v(1,t-1),mt(1,t-1));
% %     ht(1,t) = h(v(1,t-1),ht(1,t-1));
% %     nt(1,t) = n(v(1,t-1),nt(1,t-1));
% % end  
% % % plot
% % figure(104);
% % subplot(4,1,1);plot(dt:dt:tend,v(1,:),'k');hold on;
% % p_Ie = plot(dt:dt:tend,Ie+v(1,1)+10,'r');% Lowered the Ie plot to better visualization
% % ylabel('V_{mV}');
% % legend(p_Ie,sprintf('I_{Injected} = %2.2f _{\\muA}\nI_{DPP} = %2.2f _{\\muA}',I_amp,Idpp));
% % title('HH Model for Spike');
% % axis([90,160,-80,60]);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%
% % subplot(4,1,2);hold on;
% % plot(dt:dt:tend,gNat(1,:)); plot(dt:dt:tend,gKt(1,:));
% % plot(dt:dt:tend,gKt(1,:)+gNat(1,:),'k--');
% % title('Conductances');
% % legend('Na ','K','sum');
% % ylabel('G_{mS}');
% % xlim([90,160]);
% % %%%%%%%%%%%%%%%%%%%%
% % subplot(4,1,3);hold on;
% % plot(dt:dt:tend,Inat(1,:),'b');plot(dt:dt:tend,Ikt(1,:),'r');
% % plot(dt:dt:tend,Inat(1,:)+Ikt(1,:),'k--');
% % ylabel('I_{uA}');
% % legend('I_{Na}','I_K','sum');
% % title('Ionic Currents for HH model');  % Correct ionic currents are presented in page 99 Bear - Exploring the brain
% % xlim([90,160]);
% % %%%%%%%%%%%
% % subplot(4,1,4);hold on;
% % plot(dt:dt:tend,mt(1,:),'.'); plot(dt:dt:tend,ht(1,:),'--'); plot(dt:dt:tend,nt(1,:));
% % ylabel('Gating Probability');
% % legend('m- Na activation','h- Na inactivation ','n- K activataion');
% % title('Variations in Gating Variables for HH model')
% % axis([90,160,0,1]);xlabel('t_{ms}'); 
% % suptitle(sprintf('Hugkin-Huxley Model For an Active Neuron\npredicted variations for different parameters'))

%% plotting Inactivation curve for sodium channels
% % %%%%%%
% % Idpps = -4:0.1:4;   %Depolarizing Pre Pulses uA
% % 
% % v  = zeros(1,tend/dt);v_dpp  = zeros(1,tend/dt);
% % mt = zeros(1,tend/dt); ht = zeros(1,tend/dt);nt = zeros(1,tend/dt);
% % mt_dpp = zeros(1,tend/dt); ht_dpp = zeros(1,tend/dt);nt_dpp = zeros(1,tend/dt);
% % Inat = zeros(1,tend/dt); Inat_dpp = zeros(1,tend/dt);
% % 
% % %setting correct initial values
% % %calculated from previous runs
% % v(1,1:te/dt) = -64.8525*ones(1,te/dt);
% % v_dpp(1,1:te/dt) = -64.8525*ones(1,te/dt);
% % mt(1,1:te/dt) = 0.0529*ones(1,te/dt);
% % ht(1,1:te/dt) = 0.5954*ones(1,te/dt);
% % nt(1,1:te/dt) = 0.3177*ones(1,te/dt);
% % %%%
% % ina_ratio = zeros(1,length(Idpps));
% % 
% % i_count = 0;
% % 
% % for i_count = 1:length(ina_ratio)
% %     Idpp = Idpps(i_count);
% %     Ie_dpp = [zeros(1,te/dt) Idpp*ones(1,uint32((tdpp-te)/dt)) I_amp*ones(1,uint32((ts-tdpp)/dt)) zeros(1,uint32((tend-ts)/dt))];
% %     Ie = [zeros(1,te/dt) zeros(1,uint32((tdpp-te)/dt)) I_amp*ones(1,uint32((ts-tdpp)/dt)) zeros(1,uint32((tend-ts)/dt))];
% %     for t=te:tend/dt % v(0) = 0    
% %         % HH model calculations
% %         v(1,t) = V(Ie(t),v(1,t-1),mt(1,t-1),ht(1,t-1),nt(1,t-1));
% %         v_dpp(1,t) = V(Ie_dpp(t),v_dpp(1,t-1),mt_dpp(1,t-1),ht_dpp(1,t-1),nt_dpp(1,t-1));
% %         % extra plots
% %         Inat(1,t) = Ina(v(1,t-1),mt(1,t-1),ht(1,t-1));
% %         Inat_dpp(1,t) = Ina(v_dpp(1,t-1),mt_dpp(1,t-1),ht_dpp(1,t-1));
% %         %%%%%%%%%
% %         mt(1,t) = m(v(1,t-1),mt(1,t-1));
% %         ht(1,t) = h(v(1,t-1),ht(1,t-1));
% %         nt(1,t) = n(v(1,t-1),nt(1,t-1));
% %         %%%%%%%%%%%%%%
% %         mt_dpp(1,t) = m(v_dpp(1,t-1),mt_dpp(1,t-1));
% %         ht_dpp(1,t) = h(v_dpp(1,t-1),ht_dpp(1,t-1));
% %         nt_dpp(1,t) = n(v_dpp(1,t-1),nt_dpp(1,t-1));
% %     end
% %     ina_ratio(i_count) = signed_max(Inat_dpp)/signed_max(Inat);
% %     if Idpp>2
% %         ina_ratio(i_count) = -1*ina_ratio(i_count);
% %     end
% % end
% % figure(200);plot(Idpps,ina_ratio);

%%
n_frames = 1200;
F(1:n_frames) = struct('cdata',[],'colormap',[]);
fig=figure('units','normalized','outerposition',[0 0 1 1]);
ts_vector = linspace(100,149,n_frames);
sb1 = subplot(3,1,1);
sb2 = subplot(3,1,2);
sb3 = subplot(3,1,3);

I0 = 10;%uA   

myVideo = VideoWriter('HH_gatings','Uncompressed AVI');
myVideo.FrameRate = 60;  % Default 30
%myVideo.Quality = 100;    % Default 75

open(myVideo);

for j = 1:n_frames
    ts = ts_vector(j);
    Ie = [zeros(1,te/dt) I0*ones(1,uint32((ts-te)/dt)) zeros(1,uint32((tend-ts)/dt))];
    for t=te:tend/dt % v(0) = 0    
        % HH model calculations
        v(1,t) = V(Ie(t),v(1,t-1),mt(1,t-1),ht(1,t-1),nt(1,t-1));
        mt(1,t) = m(v(1,t-1),mt(1,t-1));
        ht(1,t) = h(v(1,t-1),ht(1,t-1));
        nt(1,t) = n(v(1,t-1),nt(1,t-1));
        gNat(1,t) = 1e-3*gna( mt(1,t),ht(1,t));
        gKt(1,t) = 1e-3*gk(nt(1,t));
    end  
    plot(sb1,dt:dt:tend,v(1,:));hold(sb1,'on')
    ylabel(sb1,'V_{mV}');
    legend(sb1,sprintf('I_{0} = %d \\muA',I0));
    title(sb1,'Comparsion of HH and Extended HH_{with A current} Model')
    axis(sb1,[90,160,-80,60]);
    hold(sb1,'off')
    
    plot(sb2,dt:dt:tend,gNat(1,:));hold(sb2,'on')
    plot(sb2,dt:dt:tend,gKt(1,:),'--');
    hold(sb2,'off')
    ylabel(sb2,'Conductance_{mS.cm^-2}');
    legend(sb2,'g_{Na}','g_{K}');
    title(sb2,'Variations in Na and K channel Conductance');
    axis(sb2,[90,160,0,0.04]);
    legend boxoff;
    

    plot(sb3,dt:dt:tend,mt(1,:));hold(sb3,'on')
    plot(sb3,dt:dt:tend,ht(1,:),'--');
    plot(sb3,dt:dt:tend,nt(1,:));
    hold(sb3,'off')
    ylabel('Gating Probability');xlabel('t_{ms}'); 
    legend('m- Na activation','h- Na inactivation ','n- K activataion');
    title('Variations in Gating Variables for HH model');
    axis(sb3,[90,160,0,1]);
    legend boxoff;

    %drawnow
    %F(j) = getframe(gcf);
    frame = getframe(gcf);
    writeVideo(myVideo,frame);
end
%movie(fig,F,2);
%v = VideoWriter('HH_gatings.mp4');
%open(v);
%writeVideo(F);
close(myVideo)

%movie2avi(F, 'HH_gatings.avi', 'quality', 100)