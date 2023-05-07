clear;
close all;
clc;

global g l m k

g = 9.8;
l = 10;
m = 100;
k = 150;

% vedi slide 45 - 49


syms x1 x2 real

fx = [x2; -(g/l)*sin(x1) - (k/(l*m))*x2];

Vx = m*g*l*(1-cos(x1)) + .5*l^2*m*x2^2; % Vx non è definita positiva poiché il cos(x1) vale uno per multipli di 2pigreco

% devo definire l'insieme D per cui la Vx è definita positiva,
% D = {x1 = ]-2pigreco, 2pigreco[, x2 qualsiasi}. Avendo definito D in questo modo la Vx è
% definita positiva

% ora devo vedere la derivata di Vx
Vx_punto = simplify(gradient(Vx)'*fx) % la funzione è semidefinita negativa perché x1 può essere un qualsiasi valore con x2=0

% quindi per il teorema di Lyapunov l'origine è STABILE

% ora costruiamo le circonferenze di raggio epsilon e delta

% fisso un valore di c per una curva di livello che sia appartenente a D
c = 1;
% faccio le intersezioni con gli assi così da ricavarmi le dimensioni dei semiassi
% dell'ellisse
a = max(solve(subs(Vx, x1, 0) == c, x2, 'Real', true)); % intersezione asse x1
b = max(solve(subs(Vx, x2, 0) == c, x1, 'Real', true)); % intersezione asse x2

epsilon = double(max(a,b));
delta = double(min(a,b));

% plotto i risultati ottenuti

figure
hold on
axis equal
grid on


% plot della curva di livello
[x1, x2] = meshgrid(-epsilon*1.2 : epsilon/10 : epsilon*1.2);
Vx =  m*g*l*(1-cos(x1)) + .5*l^2*m*x2.^2;
contour(x1, x2, Vx,  [c c], 'b', 'ShowText','on', 'LineWidth', 2)
% sembrano degli ellissi ma in realtà non lo sono, sembrano perché la non
% linearità è poco signficativa per i valori assunti da x1 e x2


viscircles([0,0], epsilon)
viscircles([0,0], delta)


while(true)
    x0 = ginput(1);
    scatter(x0(1), x0(2))
    grid on
    [t,y] = ode45(@non_linear, [0 40], x0);
    plot(y(:,1), y(:,2))
end

% la curva di livello blu non si vede perché epsilon e delta sono molto
% simili, basta mettere l = 1 per vedere un risultato più significativo

%%  seconda funzione di Lyapunov (slide 48)

global g l m k

g = 9.8;
l = 10;
m = 100;
k = 150;

syms x1 x2 real

fx = [x2; -(g/l)*sin(x1) - (k/(l*m))*x2];


syms p11 p12 p22 x1 x2 real
assume(p11 >= 0)
assume(p12 >= 0)
assume(p22 >= 0)

P = [p11 p12; p12 p22];
assume(det(P) > 0);

x = [x1; x2];

Vx = g*l*m*(1-cos(x1)) + .5 * x' * P * x; % funzione definita positiva
Vx_grad = (gradient(Vx)); % gradiente di VX devo prendere solo le derivate rispetto a x1 e x2 (gli ultimi due elementi del vettore)
Vx_punto = simplify(Vx_grad(4:5)' * fx); % vedi eq. slide 48

% ora devo trovare i valori di p11, p12, p22 per cui viene Vx_punto è
% definita negativa (vedi calcoli slide 48)
p12_v = k*l*.5;
p11_v = (k / (l*m)) * p12_v;
p22_v = l^2 * m;

Vx_punto = simplify(subs(Vx_punto, [p11 p12 p22], [p11_v p12_v p22_v])); % è definita negativa in D = {x1 = ]-pigreco,pigreco[ x2 qualsiasi}

% l'origine è ASINTOTICAMENTE STABILE nell'insieme D


Vx = simplify(subs(Vx, [p11 p12 p22], [p11_v p12_v p22_v])); % ora sostituisco i valori in Vx

g = 9.8;
l = 10;
m = 100;
k = 150;


c_star = subs(Vx, [x1, x2], [pi, 0]); % trovo la c*
c_star = (pi^2*k^2 + 16*g*l*m^2) / (8*m);

% ora devo prendere una curva di livello c infinitesimamente più piccola di
% c_star

c = c_star * 0.99;


% plotto i risultati ottenuti

figure
hold on
axis equal
grid on
axis([-pi*2.2 pi*2.3 -5 5])

xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})

% plot della curva di livello
[x1, x2] = meshgrid(-pi*2.5 : .5 : pi*2.5);
Vx = g*l*m*(1-cos(x1)) + .5 * (k^2/(2*m)*(x1.^2) + x1.*x2.*(k*l*.5) + x1.*x2.*(k*l*.5) + l^2*m*(x2.^2));
%contour(x1, x2, Vx,  5, 'b')%, 'ShowText','on', 'LineWidth', 2)
contour(x1, x2, Vx,  [c c], 'r', 'LineWidth', 2)
scatter([-2*pi, 0, 2*pi], [0, 0, 0], 'filled')

x0 = [2.5 .1; -2 0; 1.2 1; 1.5*pi -1; -pi 4; pi -3; -3/2*pi 1];
for i = 1 : size(x0, 1)
    scatter(x0(i, 1), x0(i, 2))
    [t, y] = ode45(@non_linear, [0 40], [x0(i, 1), x0(i, 2)]);
    plot(y(:,1), y(:,2))
end

legend('\Gamma')
% while(true)
%     x0 = ginput(1);
%     scatter(x0(1), x0(2))
%     grid on
%     [t,y] = ode45(@non_linear, [0 40], x0);
%     plot(y(:,1), y(:,2))
% end

% le curve di livello sono ruotate a causa dei termini misti

% con l=1 viene un qualcosa di più sensato

