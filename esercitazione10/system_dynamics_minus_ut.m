function dxdt = system_dynamics_minus_ut(t, x)

global k1 k2 K wn ksi

% traiettoria incognita

x1 = x(1);
x2 = x(2);
x3 = x(3);

% mi serve un filtro che mi stimi fino alla derivata seconda dato che in v
% compare la derivata del riferimento fino all'ordne due


y = x1;

r  = sin(t); % riferimento incognito

% ora la traiettoria non la conososco, quindi devo utilizzare un filtro per
% ricavarmela e ricavare anche le derivate

w1 = x(4); % stima della traiettoria
w2 = x(5); 

% stima delle derivate (implementazione del filtro vedi quad)

w1d = w2; 
w2d = -2*ksi*wn*w2 - wn^2*w1 + wn^2*r*K;

r_est = w1;
rd = w1d; % derivata della stima
rdd = w2d; % derivata seconda della stima

e = y - r_est; % errore
ed = x2 - rd; % derivata dell'errore

v = rdd - k1 * e - k2 * ed;

u = -x3^2 + v;


x1d = x2;
x2d = x3^2 + u;
x3d = -x3^3 + x1;

dxdt = [x1d; x2d; x3d; w1d; w2d];

end

