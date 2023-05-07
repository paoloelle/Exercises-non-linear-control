
clear;
close all;
clc;

global C h R E L

E = 2;
R = 3;
C = 4;
L = 5;

h = @(x) 17.76 * x - 103.79*x^2 + 229.62*x^3 + - 226.31*x^4 + 83.72*x^5;

syms x1 x2 

x1_punto = (1/C) * (-h(x1) + x2);
x2_punto = (1/L) * (-x1 - R*x2 + E); 

F = [x1_punto; x2_punto];
J = jacobian(F, [x1 x2]);

% devo trovare i punti di equilibrio del sistema
solution = solve([x1_punto == 0, x2_punto == 0], [x1 x2]);
x_eq = [vpa(solution.x1) vpa(solution.x2)]; % memorizzo i punti di equilibrio
x_eq([2, 3],:) = []; % rimuovo le righe con punti di equilibrio complessi (perchè hanno poco senso)

J_eq = zeros(2, 2, size(x_eq,1)); % preallocation for speed computation
J_eig = zeros(size(x_eq,1), 2); % preallocation for speed computation

for i = 1 : size(x_eq, 1)
    J_eq(:, : , i) = subs(J, [x1, x2], x_eq(i, :));  % jacobian evaluation in equilibrium points
    J_eig(i, :) = eig(J_eq(:, : , i)); % memorizzo gli autovalori
end

% ora guardando la tabella di pagina 45 si possono valutare i punti di
% equilibrio.
% FIXME Il codice potrebbe essere migliorato perché se cambio i valori delle
% costanti gli autovalori cambiano e non è detto che devo droppare quelle
% righe (corrispondenti ai punti di equilibrio complessi)


% plot come in slide 39
% [X, Y] = meshgrid(-1:.1:1,-1:.1:1);
% lambda1 = J_eig(1,1);
% lambda2 = J_eig(1,2);
% 
% z1 = (-1 : 0.01: 1);

% figure
% %axis equal
% xlim([-1 1])
% ylim([-1 1])
% hold on
% quiver(X,Y, lambda1*X, lambda2*Y)
% for c = [-1, .5, 1] % vedi slide 37
%     z2 = c*z1.^(lambda2/lambda1);
%     plot(z1, z2, 'LineWidth', 1)
% end



figure
hold on
s1 = scatter(x_eq([1, 3],1), x_eq([1, 3],2), 100, 'p', 'filled'); % punti di eq.stabili
s2 = scatter(x_eq(2,1), x_eq(2,2), 100, 'p', 'filled'); % punto di sella
xlim([-1 1])
ylim([-1 1])
axis equal
grid on

x0 = [.5 .1; -.3 .5; -0.9 -0.7; .5 .5; .2 .8; .7 0; .4 .7];
for i = 1 : size(x0, 1)
    scatter(x0(i, 1), x0(i, 2))
    [t, z] = ode45(@dynamics, [0 30], x0(i, :));
    plot(z(:,1), z(:,2))
end

legend([s1, s2],{'equilibrio stabile', 'sella'}, 'Location', 'southeast')


% il sistema converge ai punti di equilibrio e non converge nel punto di
% sella

