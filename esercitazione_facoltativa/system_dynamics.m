function dxdt = system_dynamics(t, x)

global K  beta ks b 


u = K * x;

u = min(abs(u), beta)*sign(u); % saturazione con soglia beta

x1d = x(2);
x2d = -ks*x(1)- b*x1d + u;

dxdt = [x1d; x2d];

end

