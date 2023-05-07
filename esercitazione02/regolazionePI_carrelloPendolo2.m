% Regolazione con informazione parziale del carrello con pendolo linearizzato

clear
clc
close all

% dati
M = 5;
m = .5;
l = 1;
Jp = m * (l^2);
g = 9.81;


A = [0      1              0               0
     0      0           -(m*g)/M           0
     0      0              0               1
     0      0       (m*g*l*(M+m))/(M*Jp)   0];

B = [0      0
     1/M    -(m*l)/(M*Jp)
     0      0
     -(m*l)/(M*Jp) (m+M)/(M*Jp)];

C = [1 0 0 0
     0 0 1 0];

D = [0 0; 0 0];



P = [0 0 0 0]'; % rumore nello stato 

omega = 1; % frequenza del riferimento
S = [0 omega; -omega 0]; % dinamica del riferimento 
S = blkdiag(S, 0);
% devo fare così perché ho un riferimento sinusoidale da inseguire ed uno
% costante quindi per il principio del modello interno devo avere la
% dinamica sia del sin che della costante nella dinamca del disturbo

Q = [0; 0]; 
% Q contiente sia il disturbo che il riferimento, solo che per il momento
% contiene solo il disturbo in uscita che è nullo, poi quando la estendo
% (con una matrice identità conterrà anche una replica del riferimento da
% inseguire)

r = size(S, 1);
n = size(A, 1);
m = size(B, 2);
q = size(C, 1);

% verifica di H1 (S antistabile)
eigen_S = eig(S); % H1 verificata

% verifica di H2 (il sistema è controllabile)
size(ctrb(A,B));
rankAB = rank(ctrb(A,B)); % H2 è verificata

% verifica di H3*
% devo costruire il sistema esteso
P = [P zeros(n, r-1)];
Ae = [A P; zeros(r, n) S];

% slide 3
Q = [-Q eye(r-1)]; % ora contiente anche la replica del riferimento
C = -C;

Ce = [C Q];
size(obsv(Ae, Ce));
rankAeCe = rank(obsv(Ae,Ce)); % H3* è verificata

% le assunzioni preliminari del teorema PI valgono, ora devo vedere se
% esistono Pigreco e Gamma. Se esistono il teorema mi dice che il problema
% PI è risolvibile

% passsaggi per ricavare Pigreco e Gamma 
L1=[-A -B; 
    C zeros(q,m)];

L2=blkdiag(eye(n), zeros(q,m));

N = [P; -Q];

Z = kron(S',L2)+kron(eye(r),L1);

sol = inv(Z)*N(:);

X=reshape(sol, [n+m,r]);

Pi=X(1:n,:);
Gamma=X(n+1:end,:);


% il problema PI è quindi risolvibile dato che ho trovato Pi e Gamma

% devo quindi ricavare le matrici F, G, H del controllore. Prima però devo
% prendere K in modo tale che il sistema A+BK abbia autovalori nel
% semipiano dx, prendo K con lqr

K = -lqr(A, B, eye(max(size(A))), 1);
eig(A+B*K); % gli autovalori sono tutti negativi del sistema a ciclo chiuso

% G devo prenderla in modo tale che il sistema Ae + G*Ce abbia gli
% autovalori nel semipiano sx
G = lqr(Ae', Ce', eye(n+r), 1);
G = G';
eig(Ae-G*Ce); % gli autovalori sono tutti negativi

G0 = G(1:n,:);
G1 = G(n+1:end, :);

F = [A-G0*C+B*K P-G0*Q+B*(Gamma - K*Pi); -G1*C S-G1*Q];

H = [K Gamma-(K*Pi)];

% ora considero i riferimenti da inseguire

Ampiezza_r = 2; % ampiezza desiderata sinusoide

syms a1 a2
% calcolo le condizioni inziali che mi servono per avere come riferimento
% quell'ampiezza (vedi slide 5)
eq = sqrt(a1^2 + a2^2) == Ampiezza_r;
solution = solve(eq, [a1 a2]);
a1 = double(solution.a1(1));
a2 = double(solution.a2(1));

Angolo_r = pi/2; % angolo desiderato

% ora lancio la simulazione
out = sim('regolazionePI_carrelloPendolo_schema.slx','StartTime', '0', 'StopTime','30');

t = out.tout; % tempo di simulazione
y1 = out.simout.Data(:,1); % posizione carrello
y2 = out.simout.Data(:,2); % angolo pendolo

figure
plot(t, [Ampiezza_r*sin(t), y1], 'LineWidth', 1);
grid on

legend('y_{1r}', 'y_1');
title('Posizione del carrello')

figure
hold on
yline(Angolo_r, 'r', 'LineWidth', 1)
plot(t,  y2, 'LineWidth', 1);
grid on
legend('y_{2r}', 'y_2');
title('Angolo del pendolo')

