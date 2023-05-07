function dxdt = non_linear(t,x)

global m k1 b1 k2 b2

dxdt = [x(2); -(k1/m)*x(1) + (k2/m) * x(1)^3 - (b1/m)*x(2) + (b2/m)*x(2)^3];

end