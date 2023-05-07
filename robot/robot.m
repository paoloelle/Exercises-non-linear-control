clear 
close all;
clc;

% dati del problema

global F1 F2 theta1 theta2 theta3 theta4 theta5 theta6 l1 l2

g = 9.81;

m1 = 5; m2 = 3;

I1 = 0.2; I2 = 0.1;

l1 = 1.5; l2 = 1;

F1 = 15; F2 = 10;

theta1 = I1 + m1*(l1^2/4) + m2*(l2^2/4 + l1^2) + I2;

theta2 = I2 + m2*(l2^2/4);

theta3 = m2 * l1 * l2/2;

theta4 = m2*l2^2/2 + I2;

theta5 = (m1*l1/2 + m2*l1)*g;

theta6 = m2*l2/2*g;

%% sistema in evoluzione libera

q0 = [0 pi/2 0 0];

[t, q] = ode45(@robot_dynamic, [0 15], q0);


G = zeros(size(t,1), 2);
for i = 1 : size(q,1)
    G(i, :) = G_q(q(i, [1,2])); % ricavo le coordinate dell'end-effector
end


% [link1_x, link1_y] = pol2cart(q(:,1), l1);
% [link2_x, link2_y] = pol2cart(q(:,1) + q(:,2), l2);

% figure
% for i = 1 : size(q,1)
%     plot([0,link1_x(i)], [0, link1_y(i)], 'LineWidth', 2)
%     grid on
%     hold on
%     plot([link1_x(i), link1_x(i) + link2_x(i)], [link1_y(i), link1_y(i) +  link2_y(i)],'LineWidth', 2, 'Color', 'g')
%     scatter( G(i, 1),  G(i, 2), 'sr', 'filled')
%     axis([-3 3 -3 3])
%     pause(.05)
%     hold off
% end

% con tau=0 e quindi senza alcun controllo sulla coppia il robot si stende
% verso il basso 

%% sistema con controllo

clc
close all


global p pd pdd kpr kdr

p = @(t) [1; sin(t)]; % coordinate desiderate end-effector
pd = @(t) [0; cos(t)]; % derivata prima coordinate desiderate end-effector
pdd = @(t) [0; -sin(t)]; % derivata seconda coordinate desiderate end-effector

kpr = diag([-10, -10]); % scelgo una matrix Hurwitz
kdr = diag([-10, -10]); % scelgo una matrix Hurwitz

global Kp Kd 

q0 = [0 pi/2 0 0 -1 pi 0 2]; % condizione iniziale 

kp1 = 10;
kp2 = 10;

Kp = diag([-kp1, -kp2]); % termine proporzinale all'errore

kd1 = 10;
kd2 = 10;

Kd = diag([-kd1; -kd2]); % termine proporzionale alla derivata dell'errore

[t, q] = ode45(@robot_dynamic_controlled, [0 20], q0);


G = zeros(size(t,1), 2);
for i = 1 : size(q,1)
    G(i, :) = G_q(q(i, [1,2]));
end


[link1_x, link1_y] = pol2cart(q(:,1), l1);
[link2_x, link2_y] = pol2cart(q(:,1) + q(:,2), l2);

% figure
% for i = 1 : size(q,1)
%     plot([0,link1_x(i)], [0, link1_y(i)], 'LineWidth', 2)
%     grid on
%     hold on
%     plot([link1_x(i), link1_x(i) + link2_x(i)], [link1_y(i), link1_y(i) +  link2_y(i)],'LineWidth', 2, 'Color', 'g')
%     scatter( G(i, 1),  G(i, 2), 'sr', 'filled')
%     axis([-3 3 -3 3])
%     pause(.05)
%     hold off
% end


figure
plot(t, [q(:, 1), q(:, 5)], 'LineWidth', 2)
grid on
legend('$q_{1}$', '$q_{1r}$', 'Interpreter', 'latex', 'FontSize', 14)

figure
plot(t, [q(:, 2), q(:,6)], 'LineWidth', 2)
grid on
legend('$q_{2}$', '$q_{2r}$', 'Interpreter', 'latex', 'FontSize', 14)


