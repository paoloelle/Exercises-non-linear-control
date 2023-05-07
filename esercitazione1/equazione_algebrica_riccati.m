% esercitazione 1: soluzione equazione di Riccati algebrica numerica

clear;
clc;
close all;

% modello massa-molla-smorzatore lineare

m = 1; % massa
c = -1; % rigidezza molla 
b = -1; % smorzamento
x0 = [1 2]'; % condizione iniziale

% il sistema diverge per valori di c e b < 0 

A = [0 1; -c/m -b/m];
B = [0; 1/m];

Q = diag([1, 1]); % costo stato
R = 1; % costo controllo

syms P11 P12 P22

P = [P11 P12; P12 P22];

ARE = transpose(A)*P + P*A - P * B * inv(R) * transpose(B)*P + Q;
solution = solve(ARE, [P11 P12 P22]);

P_manual = zeros(2,2,4); % possibile solution for ARE (le soluzioni totali sono 4 matrici 2x2)
K_manual = zeros(1, 2, 4);
eigen_manual = zeros(2, 1, 4);

% il fatto che ho 4 possibili soluzioni mi determina anche 4 possibili
% configurazioni del guadagno K e 4 possibili configurazioni di autovalori.
% Chiaramente mi prenderò la soluzione con autovalori entrambi nel
% semipiano sinistro

eigen_stability = [0; 0];

for i = 1 : 4
    
    P_manual(:, :, i) = double([solution.P11(i), solution.P12(i); solution.P12(i), solution.P22(i)]);
    K_manual(:, :, i) = inv(R) * B' * P_manual(:, :, i);
    eigen_manual(:, :, i) = eig(A - B * K_manual(:, :, i));

    if real(eigen_manual(1, :, i)) < 0 &&  real(eigen_manual(2, :, i)) < 0 % prendo gli autovalori nel semipiano sinistro
        eigen_stability = eigen_manual(:, :, i); 
        K_stability = K_manual(:, :, i); % guadagni K1, K2 che mi danno quegli autovalori
    end

end

% ora utilizzo la funzione lqr che ritorna il guadagno K, la soluzione
% dell'equazione algebrica di Riccati P, e gli autovalori (cioè tutto
% quello che ho fatto fino ad ora)

[K_lqr, P_lqr, eigen_lqr] = lqr(A, B, Q, R);

K_stability
K_lqr

eigen_stability
eigen_lqr

% sia per K che per gli autovalori ottengo lo stesso risultato tra i due
% metodi, tutto funziona correttamente.

% lqr_riccati è il modello simulink di un sistema controllato in
% retroazione con guadagno K_stability, l'obiettivo è portare entrambe le
% componenti dello stato a 0

P_lqr;
V_star = x0'*P_lqr*x0; %costo ottimo


% simulazione sistema in evoluzione libera
[t, x] = ode45(@(t, x) A*x, [0 20], x0);

figure
plot(t, [x(:, 1), x(:,2)], 'LineWidth', 1)
grid on
legend('x_1', 'x_2')
% il sistema è divergente

% simulazione sistema controllato
[t, x] = ode45(@(t, x) (A+B*(-K_lqr))*x, [0 20], x0);

figure
plot(t, [x(:, 1), x(:,2)], 'LineWidth', 1)
grid on
legend('x_1', 'x_2')




