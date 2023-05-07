function dxdt = system_dynamic(t, x)

global I1 I2 I3 k1 k2 k3


w1 = x(1);
w2 = x(2);
w3 = x(3);

w1d = 1/I1 * (I2-I3) * w3 * w2; 
w2d = 1/I2 * (I3-I1) * w3 * w1; 
w3d = 1/I3 * (I1-I2) * w1 * w2; 
        
dxdt = [w1d; w2d; w3d]; % dinamica

end