clear;
close all;
clc;

%global g l m k

g = 9.8;
l = 10;
m = 100;
k = 150;

% esercizio fatto sul quaderno 12/11 slide 77

syms x1 x2 real

fx = [x2; -(g/l)*sin(x1) - (k/(l*m))*x2];

Vx = .5 * m * l^2 * x2^2 + m * g * l * (1 - cos(x1)) % semidefinita positiva
Vx_punto = simplify(gradient(Vx)'*fx) % semidefinita negativa

% per applicare il teorema di Barbashin-Krasovski-LaSalle devo avere che la
% fx deve essere localmente Lipschitz in un dominio contenente l'origine.
% La fx rispetta questa condizione (posso utilizzare il lemma Lip.2 pag.59)
% Poi devo avere che Vx deve essere differenziabile e definita positiva,
% devo scegliere un dominio D in cui Vx è definita positiva. Il coseno vale
% 1 nei multipli di 2pigreco, posso prendere il dominio D tale per cui x1
% varia tra -2pigreco e 2pigreco esclusi. Tuttavia sapendo che il sistema
% ha punti di equilibrio instabili a multipli di pigreco prendo D come una
% curva di livello infinitesimamente più piccola di quella che interseca i
% valori di x1 a pigreco e -pigreco. Perchè poi la stima della mia regione
% di attrazione sarà contenuta in D e non posso avere una regione
% d'attrazione che include un punto di equilibrio (seppur instabile)

V_interpi = subs(Vx, [x1 x2], [pi 0]); % intersezione a pigreco
V_intermpi = subs(Vx, [x1 x2], [-pi 0]); % intersezione a meno pigreco

% Vx_punto è uguale a 0 per x2 = 0 e qualunque x1, tuttavia io ho
% considerato che x1 è compreso nell'intervallo aperto -pigreco, pigreco.
% Ora devo definire l'insieme S per e devo avere che l'unica soluzione è la
% soluzione banale x = 0. Per il principio di invarianza di LaSalle ho che
% le soluzioni convergono nel più grande invariante e devo avere che questo
% invariante è zero.
Invariant = subs(fx, x2, 0);
% quindi ottengo che nell'insieme D che ho definito l'unico invariante è
% zero, posso applicare il teorema

% ottengo quindi che l'origine è asintoticamente stabile e la curva gamma
% contenuta in D è una stima della regione di attrazione

Gamma = V_interpi * 0.99; % stima del bacino di attrazione

figure
hold on
axis equal
grid on
axis([-pi*2.2 pi*2.3 -5 5])

xticks([-2*pi -pi 0 pi 2*pi])
xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
%xline([-pi, pi], '-', {'-pi', 'pi'})


% plot della curva di livello
[x1, x2] = meshgrid(-pi*2.5 : 1 : pi*2.5);
Vx = .5 * m * l^2 * x2.^2 + m * g * l * (1 - cos(x1));

contour(x1, x2, Vx,  [Gamma Gamma], 'r', 'LineWidth', 2)
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
%     [t, y] = ode45(@non_linear, [0 40], x0);
%     plot(y(:,1), y(:,2))
% end
