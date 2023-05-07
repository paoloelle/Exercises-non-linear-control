clear;
close all;
clc;

global I1 I2 I3 k1 k2 k3

I1 = 13470;
I2 = 20450;
I3 = 27200;



%% sistema in evoluzione libera

x0 = [1 -2 -1]; % condizione iniziale

[t, y] = ode45(@system_dynamic, [0 30], x0);

figure
plot(t, [y(:,1), y(:,2), y(:,3)], 'LineWidth', 2)
grid on
legend({'$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)
%title('Dinamica del sistema in evoluzione libera')
% come da appunti il sistema è stabile infatti le traiettorie non divergono
% ma rimangono limitate in un intorno dell'origine



%% stabilizzazione con controllo lineare
% ora voglio rendere il punto di equilibrio exp. stabile. Dai calcoli si
% ottiene che basta prendere una legge di controllo in retroazione dallo
% stato

k1 = 1 * I1;
k2 = 1 * I2;
k3 = 1 * I3;

[t, y] = ode45(@dynamic_stabilization, [0 10], x0);

figure
plot(t, [y(:,1), y(:,2), y(:,3)], 'LineWidth', 2)
grid on
legend({'$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)
%title('Stabilizzazione con controllo lineare')

% dai calcoli sul quaderno ha dimosrato che l'origine è GES



%% inseguimento di traiettoria nota (feedback linearization)

% considero prima il caso in cui conosco la traiettoria e quindi anche la
% sua derivata. La traiettoria è specificata nella funzione

x0 = [1 -2 -1]; % condizione iniziale

k1 = -1;
k2 = -1;
k3 = -1;

[t, y] = ode45(@tracking_know_trajectory, [0 20], x0);

figure
%sgtitle('Inseguimento di traiettoria', 'FontSize', 11)

subplot(3,1,1)
plot(t, y(:,1), 'LineWidth', 2)
hold on
plot(t, sin(t), 'LineWidth', 2, 'LineStyle', '--')
ylim([min(x0(1), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{1}$', '$\omega_{r1}$'}, 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,1,2)
plot(t, y(:,2), 'LineWidth', 2)
hold on
plot(t, cos(t), 'LineWidth', 2, 'LineStyle', '--')
ylim([min(x0(2), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{2}$', '$\omega_{r2}$'}, 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,1,3)
plot(t, y(:,3), 'LineWidth', 2)
hold on
yline(1, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
ylim([min(x0(3), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{3}$', '$\omega_{r3}$'}, 'Interpreter', 'latex', 'FontSize', 14)



%% inseguimento di traiettoria incognita (feedback linearization)

x0 = [1 -2 -1 1 1 1]; % condizione iniziale

k1 = -1;
k2 = -1;
k3 = -1;

global alpha beta

% filter parameters
alpha = 10;
beta = alpha;

H = tf([beta], [1 alpha]); % filter fdt
figure
bode(H);
grid on


[t, y] = ode45(@tracking_unknow_trajectory, [0 20], x0);

figure
%sgtitle('Inseguimento di traiettoria incognita', 'FontSize', 11)

subplot(3,1,1)
plot(t, y(:,1), 'LineWidth', 2)
hold on
plot(t, sin(t), 'LineWidth', 2, 'LineStyle', '--')
ylim([min(x0(1), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{1}$', '$\omega_{r1}$'}, 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,1,2)
plot(t, y(:,2), 'LineWidth', 2)
hold on
plot(t, cos(t), 'LineWidth', 2, 'LineStyle', '--')
ylim([min(x0(2), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{2}$', '$\omega_{r2}$'}, 'Interpreter', 'latex', 'FontSize', 14)

subplot(3,1,3)
plot(t, y(:,3), 'LineWidth', 2)
hold on
yline(1, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
ylim([min(x0(3), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{3}$', '$\omega_{r3}$'}, 'Interpreter', 'latex', 'FontSize', 14)

% notare che ho un leggero sfasamento introdotto dal filtro. Se metti alpha
% pari ad 1 si vede come ho uno sfasamento significativo e anche il modulo
% viene attenuato coerentemente con i diagrammi di bode del filtro

%% inseguimento di traiettoria con paramtri incogniti (controllo adattativo)


% assumo di non conoscere le inerzie, quindi devo stimarle

clc

x0 = [1 -2 -1 10000 20000 30000]; % condizione iniziale


k1 =  10 * I1;
k2 =  10 * I2;
k3 =  10 * I3;

[t, y] = ode45(@tracking_unknow_inertia, [0 10], x0);


figure
subplot(3,1,1)
plot(t, y(:,1), 'LineWidth', 2)
hold on
plot(t, sin(t), 'LineWidth', 2, 'LineStyle', '--')
ylim([min(x0(1), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{1}$', '$$\omega_{r1}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location','southeast')

subplot(3,1,2)
plot(t, y(:,2), 'LineWidth', 2)
hold on
plot(t, cos(t), 'LineWidth', 2, 'LineStyle', '--')
 ylim([min(x0(2), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{2}$', '$$\omega_{r2}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location','southeast')

subplot(3,1,3)
plot(t, y(:,3), 'LineWidth', 2)
hold on
yline(1, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r');
ylim([min(x0(3), -1.5)*1.2 1.5])
grid on
legend({'$\omega_{3}$', '$$\omega_{r3}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'Location','southeast')



%% stabilizzazione sistema sottoattuato (backstepping)


% ho dimostrato nella prima sezione che l'origine del sistema è stabile,
% voglio rendere l'origine asintoticamente stabile con un sistema
% sottoattuato


% devo ricavare la formula per i valori do u1 e u2 che mi cancellano tutta
% la parte in Vpunto che mi da problemi con il segno


%-- ricavo u1 e u2 --------------------------------------------------------
clc
clear

syms I1s I2s I3s w1s w2s w3s lambda1s lambda2s lambda3s u1 u2 

w1d = 1/I1s * ((I2s-I3s) * w3s * w2s + u1);
w2d = 1/I2s * ((I3s-I1s) * w3s * w1s + u2); 
w3d = 1/I3s * (I1s-I2s) * w1s * w2s; 

w1_tilde = w1s - lambda3s*w3s;
w2_tilde = w2s - w3s^2;

w1_tilded = w1d - lambda3s*w3d;
w2_tilded = w2d - 2*w3s*w3d;


V1d = w1_tilde * w1_tilded * I1s;
V2d = w2_tilde * w2_tilded * I2s;


V1d_primo = V1d + w1_tilde * (w3s*w2_tilde + w3s^3) * (I1s - I2s);
V2d_primo = V2d + lambda3s * w2_tilde * w3s^2;

global u1_solution u2_solution

u1_solution = simplify(solve(V1d_primo == -lambda1s*w1_tilde^2, u1))
u2_solution = simplify(solve(V2d_primo == -lambda2s*w2_tilde^2, u2))

resulting_eq1 = simplify(subs(V1d_primo, u1, u1_solution)); % allora la u1_solution calcolata è corretta
resulting_eq2 = simplify(subs(V2d_primo, u2, u2_solution)); % allora la u2_solution calcolata è corretta


%------------------------------------------------------------------------


% altro procedimento per ricavare u1 e u2 (stessi risultati quindi in
% teoria vanno bene) 
% eq1 = w1_tildes*((I2s-I3s)*w2s*w3s + u1s - lambda3s*(I1s/I3s)*(I1s-I2s)*w1s*w2s) + w1_tildes*(w3s*w2_tildes + w3s^3)*(I1s-I2s);
% u1_solution = solve(eq1 == -lambda1s*w1_tildes^2, u1s);
% resulting_eq1 = simplify(subs(eq1, u1s, u1_solution)) % allora la u1_solution calcolata è corretta
% 
% 
% eq2 = w2_tildes*((I3s-I1s)*w3s*w1s + u2s - 2*w3s*(I2s/I3s)*(I1s-I2s)*w1s*w2s) + lambda3s*w2_tildes*w3s^2;
% u2_solution = solve(eq2 == -lambda2s*w2_tildes^2, u2s);
% resulting_eq2 = simplify(subs(eq2, u2s, u2_solution)) % allora la u2_solution calcolata è corretta


global lambda1 lambda2 lambda3

% lambda1, lambda2, lambda3 > 0
lambda1 = 5*10^3;
lambda2 = 5*10^3;
lambda3 = 10; 

x0 = [-2 1 3]; 

[t, y] = ode45(@system_underactuated, [0 20], x0);
figure
plot(t, [y(:,1), y(:,2), y(:,3)], 'LineWidth', 2)
grid on
legend({'$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)

% plot su un intervallo temporale molto più lungo per far vedere che anche
% omega1 va a zero
[t, y] = ode45(@system_underactuated, [0 3000], x0);
figure
plot(t, [y(:,1), y(:,2), y(:,3)], 'LineWidth', 2)
grid on
legend({'$\omega_{1}$', '$\omega_{2}$', '$\omega_{3}$'}, 'Interpreter', 'latex', 'FontSize', 14)
