function dxdt = nonlinear_dynamics(t, x)

global Cf Cr Ef Er m J lr lf Df Dr Bf Br V deltaf

vy = x(1);
r = x(2);

deltaf = 0;

alpha_f = deltaf - atan((vy + r*lf)/V);
alpha_r = -atan((vy - r*lf)/V);

f_sf = Df * sin(Cf * atan(Bf*(1-Ef)*alpha_f + Ef*atan(Bf*alpha_f)));
f_sr = Dr * sin(Cr * atan(Br*(1-Er)*alpha_r + Er*atan(Br*alpha_r)));

vyd = (f_sf)/m * cos(deltaf) + (f_sr)/m - r*V;
rd = (f_sf*lf)/J * cos(deltaf) - (f_sr*lr)/J;

dxdt = [vyd; rd];

end

