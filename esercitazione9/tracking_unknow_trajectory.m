function dxdt = tracking_unknow_trajectory(t, x)

global I1 I2 I3 k1 k2 k3 alpha beta

% desired trajectory UNKNOW (derivative too of course)
wr1 = sin(t);
wr2 = cos(t);
wr3 = 1;

w1 = x(1);
w2 = x(2);
w3 = x(3);

y1 = x(4);
y2 = x(5);
y3 = x(6);

% derivative trajectory computed with filter
y1d = -alpha * y1 + beta * wr1;
y2d = -alpha * y2 + beta * wr2;
y3d = -alpha * y3 + beta * wr3;


% errors
e1 = w1 - y1;
e2 = w2 - y2;
e3 = w3 - y3;


u1 = (-(1/I1 * (I2-I3) * w3 * w2) + y1d + k1*e1) * I1;
u2 = (-(1/I2 * (I3-I1) * w3 * w1) + y2d + k2*e2) * I2;
u3 = (-(1/I3 * (I1-I2) * w1 * w2) + y3d + k3*e3) * I3;

w1d = 1/I1 * ((I2-I3) * w3 * w2 + u1);
w2d = 1/I2 * ((I3-I1) * w3 * w1 + u2);
w3d = 1/I3 * ((I1-I2) * w1 * w2 + u3);
        
dxdt = [w1d; w2d; w3d; y1d; y2d; y3d]; % dinamica con controllo 

end