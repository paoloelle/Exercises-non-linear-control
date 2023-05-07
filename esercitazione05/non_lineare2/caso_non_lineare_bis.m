clear;
clc;
close all;

% sistema massa-molla-smorzaore caso NON lineare 2

global m k1 k2 b1 b2 

m = 100; % massa
k1 = 10; % rigidezza lineare
b1 = 50; % smorzamento lineare
k2 = 1; % rigidezza NON lineare
b2 = 5 ;% smorzamento NON lineare

syms x1 x2 real

fx = [x2; -(k1/m)*x1 + (k2/m) * x1^3 - (b1/m)*x2 - (b2/m)*x2^3];

% devo definire la funzione di Lyapunov
Vx = .5 * m * x2^2 + .5 * k1 * x1^2 - .25 * k2 * x1^4; % non è def. positiva su tutto R^2, devo trovare l'insieme D

% il plot dimostra che V(x) non è definita positiva
figure
fsurf(@(x1, x2).5 * m * x2^2 + .5 * k1 * x1^2 - .25 * k2 * x1^4, [-10 10])
xlabel('x_1')
ylabel('x_2')
zlabel('V(x)')

Vx_punto = gradient(Vx)' * fx; % -b1*x2^2 - b2*x2^4 
% Vx_punto è semidefinita negativa

% devo trovare la regione dove Vx è definita positiva, cioè l'insieme D
solution = solve(Vx == 0, [x1 x2]);
x1_Vx = double(solution.x1);
x2_Vx = double(solution.x2);

figure
scatter(x1_Vx, x2_Vx, 'filled')
grid on
axis equal

% rimuovo l'origine dall'elenco dei punti perché mi va bene che
% nell'origine si annulla
x1_Vx(x1_Vx==0) = [];
x2_Vx(x2_Vx==0) = [];

% ora devo prendere la D che non mi intercetta altri punti dove la V(x) si
% annulla

x1_min = min(abs(x1_Vx)); % prendo la c più vicina a zero 

c = double(subs(Vx, [x1, x2], [x1_min*0.99, 0])); % sarebbe il valore di Vx 

% ora trovo il valore di x2 corrispondente alla curva di livello c
x2solve(subs(Vx == c, x1, 0), x2) 

%x2_Vx_min = double(solve(subs(Vx, x1, x1_Vx_min) == c, x2))



