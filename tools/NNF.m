function Uxt = NNF(initial_Uxt,Stim,A,B,a,b,wsca,wu,wu_sign)
% Non Linear Amari Neural Field for 200 neurons
% nonlinear neural field: u' = -u + conv_x (w * I(u)) + s(x) - h
% Can change number of neurons
% A B a b parameters for mexican hat lateral inhibition intraction kernel
% wsca maxican hat kernel scaling factor
% initial_Uxt: initial Uxt should be Nx*1 row vector
% Stim: Stimulus should be 1*Nx column vector

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

Uxt = zeros(Nx,Nt);
Uxt(:,1) = initial_Uxt; % initially start with gaussian white noise

h = 1; % resting level parameter

if nargin<8, wu = 0; end


for t = 2:Nt
    for x = 1:Nx
        %computing of the convolution integral
        X = xs(x) - xs;
        Wx = wsca*((abs(X)<a)*A + ((a<=abs(X))&(abs(X)<=b))*-B);
        if wu == 1
            Wux = wu_sign*((abs(X)<b)*A)*X';
            Wx = Wx + Wux;
        end

        convWU = Wx*heaviside(Uxt(:,t-1))*dx;
        
        %computing the network
        Uxt(x,t) =  Uxt(x,t-1) + (dt/tav) * ( -Uxt(x,t-1) + convWU + Stim(x)- h );                         
    end
end

end