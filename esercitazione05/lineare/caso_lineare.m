clear;
clc;
close all;

% sistema massa-molla-smorzaore caso lineare

global m k b 

m = 100; % massa
k = 10; % rigidezza molla
b = 50; % smorzamento

syms x1 x2 real
x = [x1; x2];

fx = [0 1; -k/m -b/m] * x; % sistema linearizzato
%nonlinear = @ (x1, x2) [0 1; -k/m -b/m] * x1;%; x2];

% devo costruire la funzione di Lyapunov

P = .5 * [k 0; 0 m];
Vx = x' * P * x;
% V_x definita positiva essendo una forma quadratica con autovalori
% positivi oppure basta scriverla come .5*(m*x2^2 + kx1^2)

Vx_punto = gradient(Vx)' * fx % -50 * x2^2
% Vx_punto è semidefinita negativa perché per qualsiasi x1 mi basta che x2
% sia uguale a zero. Quindi x1 può essere anche diverso da zero 

% Quindi il sistema è STABILE

% ora devo trovare (epsilon, delta)

epsilon = 1;
beta = min(eig(P)) * epsilon^2;
delta = epsilon * sqrt(min(eig(P)) / max(eig(P)));

% plotto i risultati ottenuti

figure
hold on
axis equal
grid on

% plot ellisse
t = linspace(0,2*pi) ;
r1 = sqrt(beta/(min(eig(P)))) ; r2 = sqrt(beta/(max(eig(P)))) ;
x = r1*cos(t) ;
y = r2*sin(t) ;
plot(x,y,'r', 'LineWidth', 2)

viscircles([0,0], epsilon, 'Color', 'b')
viscircles([0,0], delta, 'Color', 'y')

text(-0.75, 0.75, '$B_\varepsilon$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(-0.35, 0, '$B_\delta$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(0.75, 0.2, '$\Omega_\beta$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')



x0 = [0.2 0.2; 0.2 -0.1; -0.15 0.1; 0 -0.2];
for i = 1 : size(x0, 1)
    scatter(x0(i, 1), x0(i, 2))
    [t, y] = ode45(@linear, [0 40], [x0(i, 1), x0(i, 2)]);
    plot(y(:,1), y(:,2))
end


% while(true)
%     x0 = ginput(1);
%     scatter(x0(1), x0(2))
%     grid on
%     [t,y] = ode45(@linear, [0 40], x0);
%     plot(y(:,1), y(:,2))
% end

% dal plot si vede che il sistema è asintoticamente stabile (potevo già
% saperlo andando a guardare gli autovalori: eig([0 1; -k/m -b/m])).


%% ora devo trovare i parametri per cui epsilon = delta

% i parametri che mi garantiscono epsilon = delta sono quelli che mi danno
% gli stessi autovalori dalla matrice P. Questo perché la matrice P è
% quella che mi definisce le dimensioni dell'ellisse di beta e quindi se
% voglio avere le due circonferenze B_delta e B_epsilon coincidenti anche
% Omega_beta deve essere una ciconferenza

clear;
clc;
close all;

% sistema massa-molla-smorzatore caso lineare

global m k b 

m = 100; % massa
k = m; % rigidezza molla
b = 50; % smorzamento

syms x1 x2
x = [x1; x2];

fx = [0 1; -k/m -b/m] * x; % sistema linearizzato
%nonlinear = @ (x1, x2) [0 1; -k/m -b/m] * x1;%; x2];

% devo costruire la funzione di Lyapunov

P = .5 * [k 0; 0 m];
Vx = x' * P * x; 
% V_x definita positiva essendo una forma quadratica con autovalori
% positivi oppure basta scriverla come .5*(m*x2^2 + kx1^2)

Vx_punto = gradient(Vx)' * fx; % -50 * x2^2
% Vx_punto è semidefinita negativa perché per qualsiasi x1 mi basta che x2
% sia uguale a zero. Quindi x1 può essere anche diverso da zero 

% Quindi il sistema è STABILE

% ora devo trovare (epsilon, delta)

epsilon = 1;
beta = min(eig(P)) * epsilon^2;
delta = epsilon * sqrt(min(eig(P)) / max(eig(P)));

% plotto i risultati ottenuti

figure
hold on
axis equal
grid on

% plot ellisse
t = linspace(0,2*pi) ;
r1 = sqrt(beta/(min(eig(P)))) ; r2 = sqrt(beta/(max(eig(P)))) ;
x = r1*cos(t) ;
y = r2*sin(t) ;
plot(x,y,'b', 'LineWidth', 2)

viscircles([0,0], epsilon)
viscircles([0,0], delta)

while(true)
    x0 = ginput(1);
    scatter(x0(1), x0(2))
    grid on
    [t,y] = ode45(@linear, [0 40], x0);
    plot(y(:,1), y(:,2))
end

% ottengo le tre circonferenze B_delta, Omega_beta e B_epsilon coincidenti
