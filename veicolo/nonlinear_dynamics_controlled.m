function dxdt = nonlinear_dynamics_controlled(t, x)

global Cf Cr Ef Er m J lr lf Df Dr Bf Br V k1 k2

vy = x(1);
r = x(2);

deltaf = -k1*vy - k2*r;

alpha_f = deltaf - atan((vy + r*lf)/V);
alpha_r = -atan((vy - r*lf)/V);

f_sf = Df * sin(Cf * atan(Bf*(1-Ef)*alpha_f + Ef*atan(Bf*alpha_f)));
f_sr = Dr * sin(Cr * atan(Br*(1-Er)*alpha_r + Er*atan(Br*alpha_r)));

vyd = (f_sf)/m * cos(deltaf) + (f_sr)/m - r*V;
rd = (f_sf*lf)/J * cos(deltaf) - (f_sr*lr)/J;

dxdt = [vyd; rd];

end


