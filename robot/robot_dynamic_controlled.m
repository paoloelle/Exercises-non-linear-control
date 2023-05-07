function dqdt = robot_dynamic_controlled(t, x)

global F1 F2 theta1 theta2 theta3 theta4 theta5 theta6 p pd pdd kpr kdr Kp Kd

% ho un sistema meccanico, considerando la singola variabile per un giunto
% avrò un termine di posizione e uno di velocità

% posizione dei giunti
q1 = x(1);
q2 = x(2);
q = [q1; q2];

% velocità dei giunti
q1d = x(3);
q2d = x(4);
qd = [q1d; q2d];

qr = [x(5); x(6)];
qrd = [x(7); x(8)]; %J_q(q)^-1 * (kpr*(G_q(q) - p(t)) + pd(t));
qrdd = J_q(q)^-1 * (pdd(t) + kdr*(J_q(q)*qrd - pd(t)) + kpr*(G_q(q) - p(t)) - diff(J_q(q))*qrd); % formula presa dal libro di sistemi robotici pag.143


% variabili di giunto (posizione angolare)


B = [theta1 + 2*theta3*cos(q2), theta2+theta3*cos(q2);
    theta2+theta3*cos(q2), theta4];

C = [-2*theta3*q2d*sin(q2), -theta3*q2d*sin(q2);
    theta3*q1d*sin(q2), 0];

g = [theta5*cos(q1)+theta6*cos(q1+q2);
    theta6*cos(q1+q2)];

F = [F1*q1d; F2*q2d];


e = q - qr;
ed = qd - qrd;
f = (C*qd + F + g); % non linearità
tau = f + B * (qrdd + Kp*e + Kd*ed); % moltiplico per B perché nell'equazione tau è pesato per B^-1
% ho un feedforward in accellerazione e due termini di feedback in
% posizione e velocità

qdd = B^(-1) * (tau - C*qd - F - g);

dqdt = [qd; qdd; qrd; qrdd];

end

