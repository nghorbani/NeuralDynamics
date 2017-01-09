function Uxt = LNF(initial_Uxt,Stim,a,b,k0)
% Linear Neural Field for 200 neurons
% LNF tav*U'(x,t) = -U(x,t) + conv(W(x),U(x,t)) + S(x,t);
% Can change number of neurons
% a b k0 parameters for Gabor function
% initial_Uxt: initial Uxt should be Nx*1 row vector
% Stim: Stimulus should be 1*Nx column vector for time
% independent or should be Nt*Nx matrix if input time dependent
% calculation of time derivative is done by forward euler formula
% calculatation of integral is done by reiman sums
% general variables setup
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

input_time_dependent = all(size(Stim)==[1,Nx]);

Uxt = zeros(Nx,Nt);
Uxt(:,1) = initial_Uxt; % initially start with gaussian white noise

% TODO: Remove these garbage variables!
c = 1; d = 2;

% simulation

for t = 2:Nt
    for x = 1:Nx
        %computing of the convolution integral
        X = xs(x) - xs;
        % CHOOSE THE INTERACTION KERNEL HERE
        % This if here is cruel!
        if input_time_dependent
            Wx = a*(exp(-X.^2/(4*b^2)).*cos(k0*X))/(sqrt(pi)*b);
        else
            Wx = exp(-c*abs(X)).*sign(X);
        end
        convWU = Wx*Uxt(:,t-1)*dx;
        
        %computing the LNF
        if input_time_dependent
            Uxt(x,t) = Uxt(x,t-1)+(dt/tav)*(-Uxt(x,t-1)+convWU+Stim(x));%1.4
        else
            Uxt(x,t) = Uxt(x,t-1)+(dt/tav)*(-Uxt(x,t-1)+convWU+Stim(x,t-1));%1.6
        end
    end
end




end