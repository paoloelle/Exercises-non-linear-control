function dqdt = robot_dynamic(t, x)

global F1 F2 theta1 theta2 theta3 theta4 theta5 theta6

% ho unn sistema meccanico, considerando la singola variabile per un giunto
% avrò un termine di posizione e uno di velocità

% variabili di giunto (posizione angolare)
q1 = x(1);
q2 = x(2);

q = [q1; q2];

% velocità dei giunti
q1d = x(3);
q2d = x(4);
qd = [q1d; q2d];


B = [theta1 + 2*theta3*cos(q2), theta2+theta3*cos(q2);
    theta2+theta3*cos(q2), theta4];

C = [-2*theta3*q2d*sin(q2), -theta3*q2d*sin(q2);
    theta3*q1d*sin(q2), 0];

g = [theta5*cos(q1)+theta6*cos(q1+q2);
    theta6*cos(q1+q2)];

F = [F1*q1d; F2*q2d];

tau = 0;

qdd = B^(-1) * (tau - C*qd - F - g);

dqdt = [qd; qdd];

end

