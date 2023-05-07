clear;
clc;
close all;

% sistema massa-molla-smorzaore caso NON lineare 1

global m k1 k2 b1 b2 

m = 100; % massa
k1 = 10; % rigidezza lineare
b1 = 50; % smorzamento lineare
k2 = 1; % rigidezza NON lineare
b2 = 5 ;% smorzamento NON lineare

syms x1 x2 

fx = [x2; -(k1/m)*x1 - (k2/m) * x1^3 - (b1/m)*x2 - (b2/m)*x2^3];

% devo definire la funzione di Lyapunov
Vx = .5 * m * x2^2 + .5 * k1 * x1^2 + .25 * k2 * x1^4; % def. positiva

Vx_punto = gradient(Vx)' * fx; % -b1*x2^2 - b2*x2^4
% Vx_punto è semidefinita negativa

% il sistema è STABILE

% ora cerco epsilon e delta, cioè le circonferenze B_epsilon e B_delta
syms c
assume(c >= 0) % dico che c è maggiore o uguale a 0 dato che la mia Vx è semidefinita positiva

% c = 5;
c = 1000; % tanto devo fissarlo un valore per c più avanti
 
% faccio le intersezioni con gli assi così da ricavarmi le dimensioni dei semiassi
% dell'ellisse
a = max(solve(subs(Vx, x1, 0) == c, x2)); % intersezione asse x1
b = max(solve(subs(Vx, x2, 0) == c, x1, 'Real', true)); % intersezione asse x2

epsilon = double(max(a,b));
delta = double(min(a,b));


% plotto i risultati ottenuti

figure
hold on
axis equal
grid on

viscircles([0,0], epsilon, 'Color', 'b')
viscircles([0,0], delta, 'Color', 'y')


% plot della curva di livello
[x1, x2] = meshgrid(-epsilon*1.2 : 0.1 : epsilon*1.2);
Vx = .5 * m * x2.^2 + .5 * k1 * x1.^2 + .25 * k2 * x1.^4;
contour(x1, x2, Vx, [c c], 'r', 'ShowText','on', 'LineWidth', 2)
% la curva di livello non è un ellisse anche se sembra

% più c è grande più la non linearità è grande e quindi più le curve di
% livello non sono ellissi


text(-0.75, 0.75, '$B_\varepsilon$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(-0.35, 0, '$B_\delta$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')

% x0 = [0.2 0.2; 0.2 -0.1; -0.15 0.1; 0 -0.2];
% for i = 1 : size(x0, 1)
%     scatter(x0(i, 1), x0(i, 2))
%     [t, y] = ode45(@non_linear, [0 40], [x0(i, 1), x0(i, 2)]);
%     plot(y(:,1), y(:,2))
% end

% while(true)
%     x0 = ginput(1);
%     scatter(x0(1), x0(2))
%     grid on
%     [t,y] = ode45(@non_linear, [0 40], x0);
%     plot(y(:,1), y(:,2))
% end

% per ottenre che le circonferenze di raggio epsilon e delta coincidono
% devo avere che a = b


