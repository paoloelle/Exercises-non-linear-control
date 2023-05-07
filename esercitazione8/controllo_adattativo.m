% controllo adattativo del sistema xd = a*x + u, a NON nota

clear;
close all;
clc;

% stabilizzazione del sistema

global a k rho

a = 2;
k = 10;
rho = 10;



%% stabilizzazione del sistema

% ho bisogon di una legge del tipo u = -a*x - k * x. Tuttavia se a non la
% conosco devo prendere una stima di a. La legge di controllo diventa u =
% -a_hat*x - k*x; dove a_hat = gamma. definisco l'errore a_tilde=a-a_hat e
% trovo il sistema xd=a_tilde-k*x; a_tilde=-gamma. Ora ricavo le funzioni
% di Lyapunov e si ottiene gamma=rho*x^2,. Applicando LaSalle si ottiene
% che l'origine è asintoticamente stabile

% x0 = [1; 1]; % condizioni iniziali
% [t, y] = ode45(@dynamic_stabilization, [0 10], x0);
% 
% x = y(:,1);
% a_hat = y(:,2);
% 
% figure
% plot(t, x, 'LineWidth', 1.5)
% grid on 
% ylim([-.1, x0(1)+.1]) % la x va a zero


figure
hold on
grid on

xlabel('$x$', Interpreter='latex')
ylabel('$\tilde{a}$', Interpreter='latex')

xlim([-4 4])
ylim([-4 4])

x0 = [-1 1; 2 -1; -2 -2; 2 -3; 1.5 2; -0.5 -2];

xline(0, '-.')
yline(0, '-.')


for i = 1 : size(x0, 1)
    
    scatter(x0(i, 1), x0(i, 2))

    [t, y] = ode45(@dynamic_stabilization, [0 10], x0(i, :));

    x = y(:,1);
    a_hat = y(:,2);

    plot(x, a_hat, 'LineWidth', 1.5)
end



%% inseguimento di traiettoria

x0 = [1; 1]; % condizioni iniziali
[t, y] = ode45(@dynamic_tracking, [0 10], x0);

x = y(:,1);

figure
plot(t, x, 'LineWidth', 1.5)
grid on

hold on
plot(t, sin(t), 'LineWidth', 1.5)
legend('x', 'x_r')





