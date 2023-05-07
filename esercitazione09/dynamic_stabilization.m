function dxdt = dynamic_stabilization(t, x)

global I1 I2 I3 k1 k2 k3


w1 = x(1);
w2 = x(2);
w3 = x(3);

u1 = -k1*w1;
u2 = -k2*w2;
u3 = -k3*w3;

w1d = 1/I1 * ((I2-I3) * w3 * w2 + u1);
w2d = 1/I2 * ((I3-I1) * w3 * w1 + u2);
w3d = 1/I3 * ((I1-I2) * w1 * w2 + u3);
        
dxdt = [w1d; w2d; w3d]; % dinamica con controllo

end