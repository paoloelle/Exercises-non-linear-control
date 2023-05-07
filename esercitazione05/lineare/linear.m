function dxdt = linear(t,x)

global m k b

dxdt = [ 0 1; -k/m -b/m] * [x(1); x(2)];

end