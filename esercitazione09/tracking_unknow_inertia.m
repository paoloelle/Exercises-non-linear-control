function dwdt = tracking_unknow_inertia(t, x)

global I1 I2 I3 k1 k2 k3

% desired trajectory
wr1 = sin(t);
wr2 = cos(t);
wr3 = 1;

% desired trajectory derivative
wr1d = cos(t);
wr2d = -sin(t);
wr3d = 0;

w1 = x(1);
w2 = x(2);
w3 = x(3);

% non conoscendo la matrix di inerzia devo utilizzare le stime
I1_hat = x(4);
I2_hat = x(5);
I3_hat = x(6);

% la legge di controllo è u = -Fw + I_hat*wrd - K*w_tilde. 

% F_wIhat = Omega * Iv_hat
Omega = [0 w2*w3 -w2*w3; -w3*w1 0 w3*w1; w1*w2 -w1*w2 0];
Iv_hat = [I1_hat; I2_hat; I3_hat];
F_wIhat = Omega * Iv_hat;

% I_hat*wrd = S * Iv_hat
S = diag([wr1d, wr2d, wr3d]);
Ihat_wrd = S * Iv_hat;

% K*w_tilde
K = diag([k1, k2, k3]);
w_tilde = [w1 - wr1; w2 - wr2; w3 - wr3];

Ihat = diag([I1_hat, I2_hat, I3_hat]);

u = -F_wIhat + Ihat_wrd - K * w_tilde; % legge di controllo

% la dinamica del sistema  è wd = (I^-1)*Fw + (I^-1)*u = (I^-1)*(Fw + u);
I = diag([I1, I2, I3]); % la I la metto per l'evoluzione del sistema anche se è incognita
Iv = [I1; I2; I3];
F_wI = Omega * Iv; % stesso discorso per la F
wd = (I^-1)*(F_wI + u); % questa è l'evoluzione del sistema


% devo far evolvere anche la stima dell'inerzia (cerco la 
% I_hatd = Gamma = (Lambda^-1) * (Omega_tilde-S)' * w_tilde
Lambda = diag([1, 1, 1]); % Lambda la scelgo arbitrariamente
w1_tilde = w1 - wr1;
w2_tilde = w2 - wr2;
w3_tilde = w3 - wr3;
Omega_tilde = [0 w2_tilde*w3_tilde -w2_tilde*w3_tilde; -w3_tilde*w1_tilde 0 w3_tilde*w1_tilde; w1_tilde*w2_tilde -w1_tilde*w2_tilde 0];
Gamma = (Lambda^-1) * (Omega_tilde-S)' * w_tilde;
Ihat_d = Gamma;


dwdt = [wd; Ihat_d];

end