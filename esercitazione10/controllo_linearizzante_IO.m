clear 
close all
clc

%% tracking di traiettoria nota (sistema con il -)

global k1 k2

k1 = 1;
k2 = 1;


x0 = [1 -2 -1];

[t, x] = ode45(@system_dynamics_minus_kt, [0 20], x0);

r = sin(t); % reference trajectory
y = x(:, 1); 

figure
plot(t, r, '--', 'LineWidth', 2);
hold on
plot(t, y, 'LineWidth', 2)
grid on
legend('r', 'y' ,'FontSize', 14)
%title('Inseguimento di traiettoria nota')

figure
plot(t, [x(:, 1), x(:, 2), x(:, 3)], 'LineWidth', 2)
grid on
legend({'$x_{1}$', '$x_{2}$', '$x_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location','southeast')
%title('Dinamica del sistema')

% la dinamica zero, ovvero la dinamica di x3 è stabile. Vedi considerazioni
% su esempio del 9/12

%% tracking di traiettoria nota (sistema con il +)


global k1 k2 



x0 = [1 -2 -1];

[t, x] = ode45(@system_dynamics_plus_kt, [0 20], x0);

r = sin(t); % reference trajectory
y = x(:, 1); 

figure
plot(t, r, '--', 'LineWidth', 2);
hold on
plot(t, y, 'LineWidth', 2)
grid on
legend('r', 'y' ,'FontSize', 14)
%title('Inseguimento di traiettoria nota')

figure
plot(t, [x(:, 1), x(:, 2), x(:, 3)], 'LineWidth', 2)
grid on
legend({'$x_{1}$', '$x_{2}$', '$x_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)
%title('Dinamica del sistema')
ylim([-3 3])
% la dinamica di x3 espolode, ma comunque sto facendo il tracking della
% traiettoria per y. ode si ferma quando x3 diverge e quindi vedo solo un
% intervallo di tempo piccolo. Se considero x3d = x3 + x1 si vede bene
% questo fenomeno perché il sistema ci mette più tempo a divergere e quindi
% vedo un intervallo di tempo più lungo in cui inseguo la traiettoria

%% tracking di traiettoria incognita (sistema con il -)


global k1 k2 K wn ksi

k1 = 1;
k2 = 1;

% parametri del filtro
K = 1;
wn = 100;
ksi = sqrt(2)/2;

% diagramma di bode del filtro
H = tf(K, [1/(wn^2) (2*ksi)/wn 1]); % filter fdt
figure
bode(H);
grid on


x0 = [1 -2 -1 1 2];

[t, x] = ode45(@system_dynamics_minus_ut, [0 20], x0);

r = sin(t); % reference trajectory
y = x(:, 1); 
r_estim = x(:, 4); % reference estimation

figure
plot(t, r, 'LineWidth', 2);
hold on
plot(t, r_estim, 'LineWidth', 2)
plot(t, y, 'LineWidth', 2)
grid on
legend('r', 'r_{e}',  'y' ,'FontSize', 11)
%title('Inseguimento di traiettoria incognita')

% se si mette wn=1 o wn=10 si vede come il filtro modifica l'ampiezza e
% fase della traiettoria stimata rispetto a quella di riferimento e come y
% insegue la traiettoria stimata dal filtro e non quella vera

figure
plot(t, [x(:, 1), x(:, 2), x(:, 3)], 'LineWidth', 2)
grid on
legend({'$x_{1}$', '$x_{2}$', '$x_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)
%title('Dinamica del sistema')
ylim([-3 3])


