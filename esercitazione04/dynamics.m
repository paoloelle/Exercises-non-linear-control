function dxdt = dynamics(t, x)

global C h R E L

x1 = x(1);
x2 = x(2);

x1d = 1/C * (-h(x1) + x2);
x2d = 1/L * (-x1 - R*x2 + E);

dxdt = [x1d; x2d];

end

