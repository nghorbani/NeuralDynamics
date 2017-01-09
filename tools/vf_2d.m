%2D vector field plotter
function vf_2d(DS,p_min,p_max,v_step,p_step)
t=0;
xx = p_min:v_step:p_max; yy = p_min:v_step:p_max;
x1=zeros(length(yy),length(xx));x2=zeros(length(yy),length(xx));
for a=1:length(xx)
    for b=1:length(yy)
        eval_f = feval(DS,t,[xx(a);yy(b)]);%evaluating derivations on each point
        x1(b,a) = eval_f(1); x2(b,a) = eval_f(2);
    end
end
norm=sqrt(x1.^2+x2.^2);
quiver(xx,yy,x1./norm,x2./norm,.5,'r');axis tight;

for x0=p_min:p_step:p_max
    for y0=p_min:p_step:p_max
        [ts,xs] = ode45(DS,[0 5],[x0 y0]);
        plot(xs(:,1),xs(:,2))
    end
end
for x0=p_min:p_step:p_max
    for y0=p_min:p_step:p_max
        [ts,xs] = ode45(DS,[0 -5],[x0 y0]);
        plot(xs(:,1),xs(:,2))
    end
end
axis([p_min p_max p_min p_max]);
