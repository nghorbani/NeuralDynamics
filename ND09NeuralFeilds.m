% Nima Ghorbani
% Goal Linear and Non-linear Neural Field

%% cleaning the workspace
clear;close all;clc;
set(0,'DefaultFigureWindowStyle','docked');

%% general values setup

tav = 10;
xend = 200; % # of neurons
tend = 100; % until time (infinity)

loc_init = -10; loc_end = 10;% initial location and infinity location

dx = (loc_end-loc_init)/xend;
dt = 0.1;

xs = loc_init+dx:dx:loc_end;
ts = dt:dt:tend;

Nx = length(xs);
Nt = length(ts);

[tt,xx] = meshgrid(ts,xs);
initial_Uxt = 0.05*randn(xend,1); %initial value (t=1) of the neural network
%initial_Uxt = zeros(xend,1);

% TODO: Remove these variables!
c = 1; d = 2;

%% LNF with only localtion dependent stimulus
% LNF tav*U'(x,t) = -U(x,t) + conv(W(x),U(x,t)) + S(x,t);

% stimulus setup
WN = normrnd(0,0.01,Nx,Nt);
Sx = exp(-xs.^2/(4*d^2))/(2*sqrt(pi)*d);
Sx_noised = Sx + WN(:,1)';

%simulation of the LNF with s(x) and s_noised(x)

%%
a = 1.0; b = 0.6; k0 = 4;
Uxt = LNF(initial_Uxt,Sx,a,b,k0);

figure(100);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('neurons');xlabel('U(x,t)');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('U(x,t)');
suptitle(sprintf('Linear Neural Field a=%2.1f, k_0=%2.1f and S_{not-noised}',a,k0));

%%
a = 1.0; b = 0.6; k0 = 4;
Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);

figure(101);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field a=%2.1f, k_0=%2.1f and S_{noised}',a,k0));

%%
a = 0.7; b = 0.6; k0 = 4;
Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);

figure(102);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field a=%2.1f, k_0=%2.1f and S_{noised}',a,k0));

%%
a = 1.5; b = 0.6; k0 = 4;
Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);

figure(103);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field a=%2.1f, k_0=%2.1f and S_{noised}',a,k0));

%%
a = 2.50; b = 0.6; k0 = 8;
Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);

figure(104);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field a=%2.1f, k_0=%2.1f and S_{noised}',a,k0));

%%
figure(105);
subplot(1,2,1); hold on;
title('U(x,\infty), k_0 = 4');
ylabel('U(x)');xlabel('neurons');
k0 = 4;b = 0.6;
for a = [0.7 1 1.5]
    Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);   
    plot(xs,Uxt(:,end));    
end
legend('a = 0.7','a = 1.0','a = 1.5')
subplot(1,2,2); hold on;
title('U(x,\infty), k_0 = 8');
ylabel('U(x)');xlabel('neurons');
k0 = 8;b = 0.6;
for a = [0.7 1 1.5]
    Uxt = LNF(initial_Uxt,Sx_noised,a,b,k0);   
    plot(xs,Uxt(:,end));    
end
legend('a = 0.7','a = 1.0','a = 1.5')

%% LNF with time and localtion dependent stimulus
% TODO a,b,k0 are not necessary here, fix for that so that there is a check

d1 = 0.5; 

% stimulus setup
v = 0.1; % traveling speed of stimulus peak
Sxt = (c/(2*sqrt(pi)*d1))*exp(-(xs'*ones(1,Nt)-ones(Nx,1)*v*ts).^2/(4*d1^2));
Uxt = LNF(initial_Uxt,Sxt); 
figure(106);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field v=%2.1f',v));


v = -0.1; % traveling speed of stimulus peak
Sxt = (c/(2*sqrt(pi)*d1))*exp(-(xs'*ones(1,Nt)-ones(Nx,1)*v*ts).^2/(4*d1^2));
%simulation of the LNF with s(x,t)
Uxt = LNF(initial_Uxt,Sxt); 
figure(107);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Linear Neural Field v=%2.1f',v));

%% 2.1
% nonlinear amari neural field: u' = -u + conv_x (w * I(u)) + s(x) - h

% stimulus setup
C = 0.6; % stimulus amplitude
d = 4; % stimulated region |x|<d
Sx = (abs(xs)<=d).*(C*(1-abs(xs)/d));
h = 1; % resting level parameter

%simulation of the Amari NNF with s(x)
%% 2.1
A = 3; B = 2; a = 1; b = 3; w_scaling = 1.0;
initial_Uxt = -h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling);

figure(108);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Non-Linear Amari Neural Field W_{scaling}=%2.1f\n W_x=-h',w_scaling));

%%% 2.1
A = 3; B = 2; a = 1; b = 3; w_scaling = 1.0;
initial_Uxt = 3.3*Sx' - h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling);

figure(109);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Non-Linear Amari Neural Field W_{scaling}=%2.1f\n W_x=-3.3*s(x)-h',w_scaling));
%%% 2.1
A = 3; B = 2; a = 1; b = 3; w_scaling = 0.1;
initial_Uxt = 3.3*Sx' - h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling);

figure(110);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle(sprintf('Non-Linear Amari Neural Field W_{scaling}=%2.1f\n W_x=-3.3*s(x)-h',w_scaling));


    
%% 2.5  wx + wu(x)
A = 3; B = 2; a = 1; b = 3; w_scaling = 1;
initial_Uxt = 3.3*Sx' - h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling,1,+1);% last two parameters wu_true and wu_sign

figure(113);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle('NNF W_x=-3.3*s(x)-h+w_u(x)');

%%% 2.5 wx - wu(x)
A = 3; B = 2; a = 1; b = 3; w_scaling = 1;
initial_Uxt = 3.3*Sx' - h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling,1,-1);% last two parameters wu_true and wu_sign

figure(114);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle('NNF W_x=-3.3*s(x)+h-w_u(x)');

%% 2.4 zero S(x)
Sx = zeros(1,Nx);

A = 3; B = 2; a = 1; b = 3; w_scaling = 1.0;
initial_Uxt = 3.3*Sx' - h * ones(Nx,1); %setting initial condition
Uxt = NNF(initial_Uxt,Sx,A,B,a,b,w_scaling);

figure(111);
subplot(3,2,[1 4]);mesh(tt,xx,Uxt);
ylabel('Neurons');xlabel('Time');zlabel('Potential');
axis('tight');
subplot(3,2,5);plot(xs,initial_Uxt);title('U(x,0)');ylabel('U(x)');xlabel('neurons');
subplot(3,2,6);plot(xs,Uxt(:,end));title('U(x,\infty)');xlabel('neurons');
suptitle('NNF with S(X)=0 and W_x=-3.3*s(x)-h');

%% 2.3 
cs = [0:0.1:0.6];
legends_temp = {};
figure(112);hold on;

for C=cs
    Sx = (abs(xs)<=d).*(C*(1-abs(xs)/d));
    initial_Uxt = 3.3*Sx' - h * ones(Nx,1);
    Uxt = NNF(initial_Uxt,Sx,A,B,a,b,1.0);
    plot(xs,Uxt(:,end));
    legends_temp = [legends_temp sprintf('C=%2.1f',C)];
end
title('NNF U(x,\infty) for W(x) =-3.3*s(x)-h with changing C in S(x)');
ylabel('U(x)');xlabel('neurons');
legend(legends_temp);

