clear;
close all;
clc;

% vedi Algoritmo di Newton 26/11

global k p p_punto l1 l2

% lunghezza dei link
l1 = 1.5;
l2 = 1;

k = diag([-1, -1]); % scelgo una matrix Hurwitz

q0 = [0; pi/2]; % angoli di partenza (condizione iniziale)
p = [(l1+l2)*cos(pi/4); (l1+l2)*sin(pi/4)]; % coordinate desiderate end-effector
p_punto = [0; 0];

[t, q] = ode45(@q_punto, [0 20], q0); %algoritmo di Newton

% traiettoria dei giunti (non richiesta)
% figure
% plot(t, q(:,1))
% 
% figure
% plot(t, q(:,2))


G = zeros(61, 2); % coord end-effector
for i = 1 : size(q)
    G(i, :) = G_q(q(i, :));
end

figure
plot(t, G(:, 1),  'LineWidth', 2) % traiettoria x end-effector
yline(p(1), 'r', 'LineWidth', 2);
grid on
legend({'$p_{1}$', '$p_{1r}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'southeast')
ylim([0 2])
%title('Traiettoria p_1 end-effector')

figure
plot(t, G(:, 2),  'LineWidth', 2) % traiettoria y end-effector
yline(p(2), 'r', 'LineWidth', 2);
grid on
legend({'$p_{2}$', '$p_{2r}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'southeast')
ylim([0 2])
%title('Traiettoria p_2 end-effector')

[link1_x, link1_y] = pol2cart(q(:,1), l1);
[link2_x, link2_y] = pol2cart(q(:,1) + q(:,2), l2);

figure

for i = 1 : size(q,1)    

    plot([0,link1_x(i)], [0, link1_y(i)], 'LineWidth', 2)
    grid on
    hold on
    plot([link1_x(i), link1_x(i) + link2_x(i)], [link1_y(i), link1_y(i) +  link2_y(i)],'LineWidth', 2, 'Color', 'g')
    scatter( G(i, 1),  G(i, 2), 'sr', 'filled')
    axis([-3 3 -3 3])
    pause
    hold off
end



%% ora gli passo una traiettoria dell'end effector
% ora avrò che p varia nel tempo

global p p_punto


% p = @(t) [0; (l1+l2)/2+(l1+l2)/4*cos(t)];
% p_punto = @(t) [0; -(l1+l2)/4*sin(t)];

q0 = [0; pi/2]; % angoli di partenza (condizione iniziale)


p = @(t) [1; sin(t)];
p_punto = @(t) [0; cos(t)];

[t, q] = ode45(@q_punto_variable, [0 20], q0); %algoritmo di Newton

% traiettoria dei giunti (non richiesta)
% figure
% plot(t, q(:,1))
% 
% figure
% plot(t, q(:,2))

G = zeros(41, 2); % coord end-effector
for i = 1 : size(q)
    G(i, :) = G_q(q(i, :));
end

figure
plot(t, G(:, 1),  'LineWidth', 2) % traiettoria x end-effector
yline(1, 'r', 'LineWidth', 2);
grid on
legend({'$p_{1}$', '$p_{1r}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'southeast')
%title('Traiettoria x end-effector')
ylim([0 2])

figure
plot(t, G(:, 2),  'LineWidth', 2) % traiettoria y end-effector
hold on
plot(t, sin(t),  'LineWidth', 2)
grid on
legend({'$p_{2}$', '$p_{2r}$'}, 'Interpreter', 'latex', 'FontSize', 14)
% title('Traiettoria y end-effector')
ylim([-2 2])

[link1_x, link1_y] = pol2cart(q(:,1), l1);
[link2_x, link2_y] = pol2cart(q(:,1) + q(:,2), l2);

figure

for i = 1 : size(q,1)    
    plot([0,link1_x(i)], [0, link1_y(i)], 'LineWidth', 2)
    grid on
    hold on
    plot([link1_x(i), link1_x(i) + link2_x(i)], [link1_y(i), link1_y(i) +  link2_y(i)],'LineWidth', 2, 'Color', 'g')
    scatter( G(i, 1),  G(i, 2), 'sr', 'filled')
    axis([-3 3 -3 3])
    pause(.05)
    hold off
end

