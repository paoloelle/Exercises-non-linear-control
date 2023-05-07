function dxdt = tracking_know_trajectory(t, x)

global I1 I2 I3 k1 k2 k3

% desired trajectory
wr1 = sin(t);
wr2 = cos(t);
wr3 = 1;

% desired trajectory derivative
wr1d = cos(t);
wr2d = -sin(t);
wr3d = 0;

w1 = x(1);
w2 = x(2);
w3 = x(3);

% errors
e1 = w1 - wr1;
e2 = w2 - wr2;
e3 = w3 - wr3;

% devo moltiplicare per le inerzie perché il controllo nelle equazioni
% della dinamica è diviso per l'inerzia
u1 = (-(1/I1 * (I2-I3) * w3 * w2) + wr1d + k1*e1) * I1;
u2 = (-(1/I2 * (I3-I1) * w3 * w1) + wr2d + k2*e2) * I2;
u3 = (-(1/I3 * (I1-I2) * w1 * w2) + wr3d + k3*e3) * I3;

w1d = 1/I1 * ((I2-I3) * w3 * w2 + u1);
w2d = 1/I2 * ((I3-I1) * w3 * w1 + u2);
w3d = 1/I3 * ((I1-I2) * w1 * w2 + u3);
        
dxdt = [w1d; w2d; w3d]; % dinamica con controllo 

end