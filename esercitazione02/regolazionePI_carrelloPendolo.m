% Regolazione con informazione parziale del carrello con pendolo linearizzato

clear
clc

% dati
M = 5;
m = .5;
l = 1;
Jp = m * l^2;
g = 9.81;

Amp = 1; % ampiezza sinusoide di riferimento
Angle = pi/2; % angolo theta di riferimento


% spazio di stato
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


% caso senza disturbo nello stato
P = zeros(length(A), 1); % matrice che "pesa" il rumore nello stato

omega = .5;
S = [omega 0; -omega 0]; % dinanimca dell'esosistema

% io devo inseguire un andamento sinusoidale per la posizione del carrello
% e un valoe costante per l'angolo del pendolo quindi S deve avere tre
% componenti, due per la sinusoide e uno per il valore costante
S = blkdiag(S, 0);

%S = 0;
Q = [-1 -1]'; %  matrice che "pesa" il rumore in uscita


%FIXME metti questi valori positivi in Q e cambia il segno della somma in simulink

%TODO qui ho messo una dinamica del disturbo nulla, quindi ho un disturbo
%costante e quindi devo avere che inseguno un riferimento costante da
%inseguire (principio del modello interno). Oppure metti una matrice S che
%definisce una dinamica sinusoidale. IL problema vuole che i

r = size(S, 1);
n = size(A, 1);
m = size(B, 2);
q = size(C, 1);



% devo verifcare che valgono le assunzioni H1, H2, H3* per provare ad
% applicare il teorema PI

% verifica di H1
lambda_S = eig(S); % spettro di S
% ho un solo autovalore in 0 e quindi vale H1


% verifica di H2
ctrb_AB = ctrb(A,B); % matrix di controllabilità
rankCtrb_AB = rank(ctrb_AB); % il rango è massimo quinidi il sistema è controllabile e quindi vale H2


% verifica di H3*
% costruisco il sistema esteso

P = [P zeros(n, r-1)]

Ae = [A P;
    zeros(r, n) S]

Q = [Q eye(r-1)];

Ce = [C Q]

obsv_AeCe = obsv(Ae, Ce); % matrix di osservabilità costrutita da Ae, Ce 
rankObsv_AeCe = rank(obsv_AeCe)
% la matrice di osservabilità ha rango pieno, quindi il sistema (Ae, Ce) è
% osservabile e quindi anche rilevabile. Quindi vale l'assuzione H3*

% le assunzioni preliminari del teorema PI valgono, ora devo vedere se
% esistono Pigreco e Gamma. Se esistono il teorema mi dice che il problema
% PI è risolvibile

% ora devo trovare le matrici Pigreco e Gamma

syms pi1 pi2 pi3 pi4 g1 g2
Pigreco = [pi1; pi2; pi3; pi4];
Gamma = [g1; g2];

eq2 = C*Pigreco + Q;
solution2 = solve(eq2, Pigreco);
Pigreco = double([solution2.pi1; solution2.pi2; solution2.pi3; solution2.pi4])
eq1 = A*Pigreco + B*Gamma + P - Pigreco*S;
solution1 = solve(eq1, Gamma);
Gamma = double([solution1.g1; solution1.g2])

% esistono le matrici Pigreco e Gamma, quindi il problema PI è risolubile
% per il teorema PI

% dal teorema CI posso ricavare L = Gamma - K*Pigreco. K la prenderò in
% modo tale che gli autovalori di A+BK siano nel semipiano sinistro

K = -lqr(A, B, eye(max(size(A))), 1)
eig(A+B*K) % gli autovalori sono tutti negativi del sistema a ciclo chiuso

L = Gamma - K * Pigreco

% ora scrivo la legge di controllo

H = [K L] 

% ora costruisco G
G = lqr(Ae', Ce', eye(n+r), 1);
G = G';

eig(Ae-G*Ce)

G0 = G(1:n,:);
G1 = G(n+1:end, :);


F = [A-G0*C+B*K P-G0*Q+B*(L); -G1*C S-G1*Q];




