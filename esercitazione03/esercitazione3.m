%% Soluzione numerica di equazioni differenziali non lineari 
 
clc 
clear 
close all  

% in riferimento ai teoremi delle slide io avrò sempre che t0 = 0 e t sarà
% sempre maggiore di zero per avere un senso, cioè ho il tempo che aumenta

%% equazione 1

x0 = -1; % condizione iniziale
[t, x] = ode45(@(t,x) x^2, [0 20], x0); 

figure
plot(t, x, 'LineWidth', 1);
xlabel('t'); 
ylabel('x(t)');
grid on

[t, x] = ode45(@(t,x) x^2, [0 20], 1);
hold on 
plot(t, x, 'LineWidth', 1);
xlabel('t'); 
ylabel('x(t)');

%xlim([-5 5])
ylim([-2 5])
legend('x_0 = -1', 'x_0 = 1')

%% equazione 2

x0 = 1; % condizione iniziale
[t, x] = ode45(@(t,x) -x^3, [0 20], x0); 

figure
plot(t, x, 'LineWidth', 1);
xlabel('t'); 
ylabel('x(t)');
grid on

%% equazione 3

x0 = 1; % condizione iniziale
[t, x] = ode45(@(t,x) x^3, [0 20], x0);
[t2, x2] = ode45(@(t,x) x^3, [0 -20], x0); 
 
figure
hold on
plot(t, x, 'LineWidth', 1);
plot(t2, x2, 'LineWidth', 1);
scatter(t(1),x(1), 'filled', 'og')
xlabel('t'); 
ylabel('x(t)');
grid on
ylim([0 10])
legend('t_f = 20', 't_f = -20', 'x_0',  'Location','northwest')