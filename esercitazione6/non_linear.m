function dxdt = non_linear(t,x)

global g l m k

dxdt = [x(2); -(g/l)*sin(x(1)) - (k/(l*m))*x(2)];

end