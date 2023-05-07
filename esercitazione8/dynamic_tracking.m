function dxdt = dynamic_tracking(t, x0)

global a k rho

x = x0(1); 
a_hat = x0(2); 

xr = sin(t); % traiettoria di riferimento
xr_d = cos(t);

e = x - xr;

u = -(a_hat*x) + xr_d - k*e;
xd = a*x + u;

a_hat_d = rho* e*x;

dxdt = [xd; a_hat_d];

end

