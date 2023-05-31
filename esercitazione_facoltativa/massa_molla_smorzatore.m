clear;
close all;
clc;


yalmip('clear')

m = 10; % massa
% ks rigidezza (costante elastica della molla) > 0
% b  smorzamento (attrito viscoso) > 0

% ho il sistema massa molla smorzatore nel caso lineare

syms ks b real % ks e b li prendo positivi altrimenti non avrebbero senso

A = [0 1; -ks/m -b/m];

%% analisi di stabiità esatta

% devo vedere le proprietà di stabilità al variare di ks e b
% ora non ho il controllo dato che sto facendo analisi della stabilità del
% sistema.

syms ks b positive % ks e b li prendo positivi altrimenti non avrebbero senso

A = [0 1; -ks/m -b/m];
% ho il sistema lineare, devo scegliere un qualche modo per ottenere le
% proprietà di stabilità. Uso il criterio di Routh.
pc = charpoly(A); % coefficenti del polinomio caratteristico

% Per il criterio di Routh i coefficenti del polinomio caratteristico
% devono avere tutti lo stesso segno, quindi devo avere che b e ks devono
% essere maggiori di 0 per avere la stabilità asintotica e infatti ha senso
% poichè ks e b sono termini dissipativi e nelle equazioni della dinamica
% del sistema hanno un meno davanti. Vedi teorema 3.6 pag.68 libro
% fondamenti di controlli automatici

%% analisi di stabilità robusta (incertezza strutturata)

% ora considero di avere un'incertezza sui parametri ks e b. Considero il
% caso di incertezza politopica. Quindi so che ks e b variano in un range
% di valori. Anche qui non ho il controllo.

% A'T + PA < 0 and P > 0 devo trovare una P che soddisfa queste condizioni
% per avere stabilità asintotica (Lyapunov). Quindi devo ottenere che tutti
% gli autovalori di G sono negativi per ogni possibile valore di ks e di b
% nell'intervallo di incertezza (teorema per la stabilità quadratica).

n = size(A,1);
P = sdpvar(n,n);

ks_min = 1;
ks_max = 10;

b_min = 1;
b_max = 10;

% considero i vertci dell'insieme determinati dall'intervallo di variazione
% dell'incertezza
A1 = double(subs(A, [ks b], [ks_min b_min]));
A2 = double(subs(A, [ks b], [ks_min b_max]));
A3 = double(subs(A, [ks b], [ks_max b_min]));
A4 = double(subs(A, [ks b], [ks_max b_max]));


G = [P>=eye(n), (A1'*P)+(P*A1) <= -eye(n)
     P>=eye(n), (A2'*P)+(P*A2) <= -eye(n)
     P>=eye(n), (A3'*P)+(P*A3) <= -eye(n)
     P>=eye(n), (A4'*P)+(P*A4) <= -eye(n)];

options=sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
solution = solvesdp(G,0,options);

P = value(P); 

% trovo la P che soddisfa le disuguaglianze per i vertici dell'insieme
% dell'incertezza e quindi il sistema è quadraticamente stabile


%% sintesi del controllore

% vedi esempio 9/12 

B = [0; 1/m]; % il controllo agisce sull'accelerazione

% adesso verifico se il sistema è stabilizzabile, se così fosse esiste una
% matrix P che poi mi consente di ricavare il guadagno K del controllore

% costruisco la matrice di controllabilità
ctrAB = [B A*B];
det(ctrAB)
% la matrix ha il determinante diverso da zero e quindi ha rango pieno e
% quindi la coppia (A,B) è controllabile (e quindi stabilizzabile)
% indipendentemente dai valori di ks e b

% ora devo definire i vertici dell'insieme convesso basandomi sui limiti
% dell'incertezza di ks e b

ks_min = 1;
ks_max = 10;

b_min = 1;
b_max = 10;

A1 = double(subs(A, [ks b], [ks_min b_min]));
A2 = double(subs(A, [ks b], [ks_min b_max]));
A3 = double(subs(A, [ks b], [ks_max b_min]));
A4 = double(subs(A, [ks b], [ks_max b_max]));

% definisco la LMI

n = size(A,1);
m = size(B,2);

W = sdpvar(n,n, 'symmetric');

Y = sdpvar(m,n, 'full');

% considero i vertici dell'insieme
F = [W>=eye(n), (A1*W+B*Y)+(A1*W+B*Y)'<= -eye(n)
     W>=eye(n), (A2*W+B*Y)+(A2*W+B*Y)'<= -eye(n)
     W>=eye(n), (A3*W+B*Y)+(A3*W+B*Y)'<= -eye(n)
     W>=eye(n), (A4*W+B*Y)+(A4*W+B*Y)'<= -eye(n)];

options=sdpsettings('solver', 'sedumi', 'sedumi.eps', 1e-8, ...
                'sedumi.cg.qprec', 1, 'sedumi.cg.maxiter', 49, ...
                'sedumi.stepdif', 2);
solution = solvesdp(F,0,options);

P0 = inv(value(W));
K = value(Y) * inv(value(W)); % K stabilizzante



%% simulazione del sistema con i guadagni K trovati e con saturazione

global K beta ks b 


% beta è il parametro per la funzione saturazione


x0 = [2 -3];

ks = 1;
b = 1;

beta = 10;
[t1, x1] = ode45(@system_dynamics, [0 10], x0);

figure
plot(t1, [x1(:, 1), x1(:,2)], 'LineWidth',1)
grid on
legend('x_1', 'x_2')
%title('\beta = 10')

beta = .1;
[t2, x2] = ode45(@system_dynamics, [0 10], x0);

figure
plot(t2, [x2(:, 1), x2(:,2)], 'LineWidth',1)
grid on
legend('x_1', 'x_2')
%title('\beta = 0.1')

A = [0 1; -ks/m -b/m];
B = [0; 1/m];

eig(A+B*K);

%% plot della funzione saturazione con andamento del controllo

u1 = K*x1';
beta = 10;
u1 = min(abs(u1), beta) .* sign(u1);

figure
hold on
plot(t1, u1, 'LineWidth', 1)
yline(beta, 'LineWidth', 1, 'Color', 'r', 'LineStyle','--')
yline(-beta, 'LineWidth', 1, 'Color', 'r', 'LineStyle','--')
grid on
ylim([-15 15])
legend('u', '\beta')

u2 = K*x2';
beta = .1;
u2 = min(abs(u2), beta) .* sign(u2);

figure
hold on
plot(t2, u2, 'LineWidth', 1)
yline(beta, 'LineWidth', 1, 'Color', 'r', 'LineStyle','--')
yline(-beta, 'LineWidth', 1, 'Color', 'r', 'LineStyle','--')
grid on
ylim([-.15 .15])
legend('u', '\beta')
