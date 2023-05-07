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
Vx = .5 * m * x2^2 + .5 * k1 * x1^2 - .25 * k2 * x1^4; % non è def. positiva su tutto R^2

Vx_punto = gradient(Vx)' * fx; % -b1*x2^2 - b2*x2^4 
% Vx_punto è semidefinita negativa

% devo trovare la regione dove Vx è definita positiva, cioè l'insieme D

% trovo gli equilibri del sistema
x_eq = solve([fx(1) == 0, fx(2) == 0], [x1, x2]);
x1_eq = (x_eq.x1)
x2_eq = (x_eq.x2)

% ora prendo una curva di livello che passa per gli equilibri trovati
c = subs(Vx, [x1, x2], [x1_eq(2) , x2_eq(2)]);

x1_inter = [x1_eq(2); x1_eq(3)]; % intersezioni con l'asse x1 (la uso per i plot)
% mi basta fare così per trovare l'intersezione con l'asse x1 perché x2 è
% uguale a zero per i punti x1_eq(2) e x1_eq(3)

% cerco l'intersezione con l'asse x2 per la curva di livello Vx trovata che
% intercetta i punti di equilibrio trovati diversi dall'origine
Vx_intx2 = subs(Vx, x1, 0);
x2_inter = solve(Vx_intx2 == c, x2) % intersezioni con l'asse x2 (usata per i plot)

x1_min = double(min(x1_inter));
x1_max = double(max(x1_inter));

x2_min = double(min(x2_inter));
x2_max = double(max(x2_inter));

figure
fsurf(@(x1, x2) .5 * m * x2^2 + .5 * k1 * x1^2 - .25 * k2 * x1^4, [x1_min x1_max x2_min x2_max])
xlabel('x_1')
ylabel('x_2')
zlabel('V(x)')


% ho trovato l'insieme D.
% In questo insieme D ho la stabilità dell'origine,
% ora devo fare in modo che omega_c e D abbiano
% intersezione nulla (slide 57-58, in queste slide parla di bacino di
% attrazione perché ha la stabiltà asintotica)

c = c - 0.001; % sarebbe la omega_c

a = max(solve(subs(Vx, x1, 0) == c, x2, 'Real', true)); % intersezione asse x1
b = max(solve(subs(Vx, x2, 0) == c, x1, 'Real', true)); % intersezione asse x2


epsilon = double(max(a,b));
delta = double(min(a,b));

figure
hold on
axis equal
grid on

% plot della curva di livello
[x1, x2] = meshgrid(-epsilon*1.2 : 0.1 : epsilon*1.2);
Vx = .5 * m * x2.^2 + .5 * k1 * x1.^2 - .25 * k2 * x1.^4;
contour(x1, x2, Vx,  [c c], 'r', 'LineWidth', 2);
% la curva di livello non è un ellisse anche se sembra

viscircles([0,0], epsilon, 'Color', 'b');
viscircles([0,0], delta, 'Color', 'y');

% plot dell'insieme D
x = [x1_inter(1), x1_inter(2), x1_inter(2), x1_inter(1), x1_inter(1)];
y = [x2_inter(1), x2_inter(1), x2_inter(2), x2_inter(2), x2_inter(1)];
plot(x, y, 'g', 'LineWidth', 2);

text(-0.75, 0.75, '$B_\varepsilon$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(-0.35, 0, '$B_\delta$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(-1, 0.75, '$\Omega_\beta$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')
text(1, 0, '$D$', 'BackgroundColor','w',  'FontSize', 15, 'Interpreter', 'latex')

x0 = [.5 .2; 0.1 .6; -0.4 0.2; 0 -0.45];
for i = 1 : size(x0, 1)
    scatter(x0(i, 1), x0(i, 2))
    [t, y] = ode45(@non_linear, [0 40], [x0(i, 1), x0(i, 2)]);
    plot(y(:,1), y(:,2))
end

% while(true)
%     x0 = ginput(1);
%     scatter(x0(1), x0(2))
%     grid on
%     [t,y] = ode45(@non_linear, [0 40], x0);
%     plot(y(:,1), y(:,2))
% end


% per ottenre che le circonferenze di raggio epsilon e delta coincidono
% devo avere che a = b






