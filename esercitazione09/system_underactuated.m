function dxdt = system_underactuated(t, x)


w1 = x(1);
w2 = x(2);
w3 = x(3);

global I1 I2 I3 lambda1 lambda2 lambda3 


u1 = -(I3*lambda1*w1 - I3^2*w2*w3 - I1^2*lambda3*w1*w2 + I1*I3*w2*w3 - I3*lambda1*lambda3*w3 + I1*I2*lambda3*w1*w2)/I3;
u2 = -(lambda2*(- w3^2 + w2)^2 - I2*(- w3^2 + w2)*((w1*w3*(I1 - I3))/I2 + (2*w1*w2*w3*(I1 - I2))/I3) + lambda3*w3^2*(- w3^2 + w2))/(- w3^2 + w2);


w1d = 1/I1 * ((I2-I3) * w3 * w2 + u1);
w2d = 1/I2 * ((I3-I1) * w3 * w1 + u2); 
w3d = 1/I3 * (I1-I2) * w1 * w2;


dxdt = [w1d; w2d; w3d];


end

