function dxdt = system_dynamics_minus_kt(t, x)

global k1 k2

% traiettoria nota

x1 = x(1);
x2 = x(2);
x3 = x(3);

y = x1;


r  = sin(t); % riferimento

rd = cos(t); %derivata prima di r rispetto al tempo
rdd = -sin(t); %derivata seconda di r rispetto al tempo

e = y - r; % errore
ed = x2 - rd; % derivata dell'errore

v = rdd - k1 * e - k2 * ed;

u = -x3^2 + v;


x1d = x2;
x2d = x3^2 + u;
x3d = -x3^3 + x1;

dxdt = [x1d; x2d; x3d];


end

