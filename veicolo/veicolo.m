clear 
close all
clc

global Cf Cr Ef Er m J cr cf lr lf l g Df Dr Bf Br V

Cf = 1.3;
Cr = 1.3;
Ef = -1.9990;
Er = -1.7908;
m = 1864;
J = 3654;
cr = 213800;
cf = 101600;
lr = 1.32;
lf = 1.51;
l = lr + lf;
g = 9.81;
Df = (1.2 * m * g * lr)/l;
Dr = (1.2 * m * g * lf)/l;
Bf = cf / (Cf*Df);
Br = cr / (Cr*Dr);



%% analisi della regione di attrazione del sistema al variare della velocità V

% considero tre velocità in m/s V1 = 15 (54 km/h)
% V2 = 25 (90 km/h), V3 = 35 (126 km/h) V4 = 50 (180km/h)
% l'evoluzione del sistema dipende dalla velocità

V1 = 15;
V2 = 25;
V3 = 35;
V4 = 50;

x0 = [-2.5 .7; -1.5 -1.5; 2 2; 3.5 3.5]; % condizioni inizali


figure
hold on
grid on
axis([-6 6 -5 5])

scatter(x0(:,1), x0(:,2), 'black')

for i = 1 :size(x0,1)
       
    V = V1;
    [t, y] = ode45(@nonlinear_dynamics, [0 40], x0(i,:));
    plot(y(:,1), y(:,2), 'r', 'LineWidth', 1)
     

    V = V2;
    [t, y] = ode45(@nonlinear_dynamics, [0 40], x0(i,:));
    plot(y(:,1), y(:,2), 'b', 'LineWidth', 1)
       

    V = V3;
    [t, y] = ode45(@nonlinear_dynamics, [0 40], x0(i,:));
    plot(y(:,1), y(:,2), 'y', 'LineWidth', 1)

    V = V4;
    [t,y] = ode45(@nonlinear_dynamics, [0 40], x0(i,:));
    plot(y(:,1), y(:,2), 'g', 'LineWidth', 1)
    
end

legend('Condizione inziale', '$V_{1}=54 \frac{km}{h}$', '$V_{2}=90 \frac{km}{h}$','$V_{3}=126 \frac{km}{h}$', '$V_{4}=180 \frac{km}{h}$', 'Interpreter', 'latex', 'Location','southeast')
%title('Dinamica del sistema in evoluzione libera')

% all'aumentare della velocità diminuisce il bacino di attrazione


%% Analisi della linearizzazione del sistema nell'intorno dell'origine

% devo prima linearizzare il sistema

syms vys rs Vs real

deltaf = 0;

alpha_f = deltaf - atan((vys + rs*lf)/Vs);
alpha_r = -atan((vys - rs*lf)/Vs);

f_sf = Df * sin(Cf * atan(Bf*(1-Ef)*alpha_f + Ef*atan(Bf*alpha_f)));
f_sr = Dr * sin(Cr * atan(Br*(1-Er)*alpha_r + Er*atan(Br*alpha_r)));

vyd = (f_sf)/m * cos(deltaf) + (f_sr)/m - rs*Vs;
rd = (f_sf*lf)/J * cos(deltaf) - (f_sr*lr)/J;


A = jacobian([vyd,rd], [vys, rs]); % compute linearization
A = subs(A, [vys, rs], [0, 0]) % evaluate in (0,0)

% ho il sistema linearizzato nell'intorno dell'origine

% considero tre velocità in m/s V1 = 15 (54 km/h)
% V2 = 25 (90 km/h), V3 = 35 (126 km/h)

V1 = 15;
A1 = subs(A, Vs, V1);
e1 = eig(A1);


V2 = 25;
A2 = subs(A, Vs, V2);
e2 = eig(A2);

V3 = 35;
A3 = subs(A, Vs, V3);
e3 = eig(A3);

V4 = 50;
A4 = subs(A, Vs, V4);
e4 = eig(A4);

figure

hold on
plot(real(e1),imag(e1),'*r', 'LineWidth', 1)
plot(real(e2),imag(e2),'*b', 'LineWidth', 1)
plot(real(e3),imag(e3),'*y', 'LineWidth', 1)
plot(real(e4),imag(e4),'*g', 'LineWidth', 1)

xline(0, ':')
yline(0, ':')

xlabel('Real')
ylabel('Imaginary')
legend('$V_{1}=54 \frac{km}{h}$', '$V_{2}=90 \frac{km}{h}$','$V_{3}=126 \frac{km}{h}$', '$V_{4}=180 \frac{km}{h}$', 'Interpreter', 'latex')
xlim([-inf, 1])
ylim([-10, 10])
axis equal
grid on
%title('Autovalori del sistema in evoluzione libera')

% per le quattro velocità scelte gli autovalori sono a parte reale
% negativa, quindi per tutte e quattro le velocità l'origine è localmente
% asintoticamente stabile per il sistema non lineare in realtà per
% l'estensione del teorema di linearizzazione di Lyapunov è anche
% localmente esponenzialmente stabile


%% tuning legge di controllo  

% prima considero il sistema lineare controllato e mostro che anche con k1
% e k2 il sistema non lineare è esponenzialmente stabile (e quindi anche
% asintoticamente) poi vedo come si comporta il sistema non lineare
% controllato con i guadagni di k1 e k2 che ho trovato cioè vedi il bacino
% di attrazione se cambia prendendo le condizioni iniziali che hai usato
% prima


% devo prima linearizzare il sistema controllato

syms vys rs Vs real
syms k1s k2s

deltaf = -k1s*vys - k2s*rs;

alpha_f = deltaf - atan((vys + rs*lf)/Vs);
alpha_r = -atan((vys - rs*lf)/Vs);

f_sf = Df * sin(Cf * atan(Bf*(1-Ef)*alpha_f + Ef*atan(Bf*alpha_f)));
f_sr = Dr * sin(Cr * atan(Br*(1-Er)*alpha_r + Er*atan(Br*alpha_r)));

vyd = (f_sf)/m * cos(deltaf) + (f_sr)/m - rs*Vs;
rd = (f_sf*lf)/J * cos(deltaf) - (f_sr*lr)/J;


A = jacobian([vyd,rd], [vys, rs]); % compute linearization
A = subs(A, [vys, rs], [0, 0]); % evaluate in (0,0)

V1 = 15;
A1 = subs(A, Vs, V1);
e1 = eig(A1);

V2 = 25;
A2 = subs(A, Vs, V2);
e2 = eig(A2);

V3 = 35;
A3 = subs(A, Vs, V3);
e3 = eig(A3);

V4 = 50;
A4 = subs(A, Vs, V4);
e4 = eig(A4);

% gli autovalori che ho trovato per le velocità dipendono da k1s e k2s

solution1 = solve(e1 == [-20; -20], [k1s, k2s]);
A1c = subs(A1, [k1s, k2s], [solution1.k1s, solution1.k2s]);

solution2 = solve(e2 == [-20; -20], [k1s, k2s]);
A2c = subs(A2, [k1s, k2s], [solution2.k1s, solution2.k2s]);

solution3 = solve(e3 == [-20; -20], [k1s, k2s]);
A3c = subs(A3, [k1s, k2s], [solution3.k1s, solution3.k2s]);

solution4 = solve(e4 == [-20; -20], [k1s, k2s]);
A4c = subs(A4, [k1s, k2s], [solution4.k1s, solution4.k2s]);


figure
hold on
grid on


plot(real(eig(A1c)), imag(eig(A1c)),'*r', 'LineWidth', 1)
plot(real(eig(A2c)), imag(eig(A2c)),'*b', 'LineWidth', 1)
plot(real(eig(A3c)), imag(eig(A3c)),'*y', 'LineWidth', 1)
plot(real(eig(A4c)), imag(eig(A4c)),'*g', 'LineWidth', 1)

xline(0, ':')
yline(0, ':')

xlabel('Real')
ylabel('Imaginary')
legend('$V_{1}=54 \frac{km}{h}$', '$V_{2}=90 \frac{km}{h}$','$V_{3}=126 \frac{km}{h}$', '$V_{4}=180 \frac{km}{h}$', 'Interpreter', 'latex')
xlim([-inf, 1])
ylim([-10, 10])
axis equal
grid on
%title('Autovalori del sistema controllato')

% ho tutti gli autovalori in -20, i gudagni cambieranno a seconda della
% velocità per portare gli autovalori in -20

% plot dei guadagni
K = [double(solution1.k1s) double(solution1.k2s); 
     double(solution2.k1s) double(solution2.k2s);
     double(solution3.k1s) double(solution3.k2s);
     double(solution4.k1s) double(solution4.k2s);];

figure
hold on
grid on

x = categorical({'V_1','V_2','V_3','V_4'});


b1 = bar(x(1), K(1,:), 'FaceColor','r');
b2 = bar(x(2), K(2,:), 'FaceColor','b');
b3 = bar(x(3), K(3,:), 'FaceColor','y');
b4 = bar(x(4), K(4,:), 'FaceColor','g');


text(b1(1).XEndPoints, b1(1).YEndPoints, 'k1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(b1(2).XEndPoints, b1(2).YEndPoints, 'k2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

text(b2(1).XEndPoints, b2(1).YEndPoints, 'k1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(b2(2).XEndPoints, b2(2).YEndPoints, 'k2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

text(b3(1).XEndPoints, b3(1).YEndPoints, 'k1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(b3(2).XEndPoints, b3(2).YEndPoints, 'k2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

text(b4(1).XEndPoints, b4(1).YEndPoints, 'k1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(b4(2).XEndPoints, b4(2).YEndPoints, 'k2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

legend('$V_{1}=54 \frac{km}{h}$', '', '$V_{2}=90 \frac{km}{h}$', '', '$V_{3}=126 \frac{km}{h}$',  '', '$V_{4}=180 \frac{km}{h}$',  '', ...
    'Interpreter', 'latex', 'Location','northwest')


%title('Guadagni del sistema controllato')


% prima avevo degli autovalori complessi, quindi non è facile confrontarli
% con degli autovalori reali e per questo motivo che il guadagno k2 ci sta
% che venga negativo (ricorda che deltaf ha già il segno meno)

global k1 k2

x0 = [-2.5 .7; -1.5 -1.5; 2 2; 3.5 3.5]; % condizioni inizali

figure
hold on
grid on

scatter(x0(:,1), x0(:,2), 'black')

for i = 1 :size(x0,1)
       
    k1 = double(solution1.k1s);
    k2 = double(solution1.k2s);
    [t, y1] = ode45(@nonlinear_dynamics_controlled, [0 40], x0(i,:));
    plot(y1(:,1), y1(:, 2), 'r', 'LineWidth', 1)
     

    k1 = double(solution2.k1s);
    k2 = double(solution2.k2s);
    [t, y2] = ode45(@nonlinear_dynamics_controlled, [0 40], x0(i,:));
    plot(y2(:,1), y2(:,2), 'b', 'LineWidth', 1)
       

    k1 = double(solution3.k1s);
    k2 = double(solution3.k2s);
    [t, y3] = ode45(@nonlinear_dynamics_controlled, [0 40], x0(i,:));
    plot(y3(:,1), y3(:,2), 'y', 'LineWidth', 1)

    k1 = double(solution4.k1s);
    k2 = double(solution4.k2s);
    [t,y4] = ode45(@nonlinear_dynamics_controlled, [0 40], x0(i,:));
    plot(y4(:,1), y4(:,2), 'g', 'LineWidth', 1)
    
end

legend('Condizione inziale', '$V_{1}=54 \frac{km}{h}$', '$V_{2}=90 \frac{km}{h}$','$V_{3}=126 \frac{km}{h}$', '$V_{4}=180 \frac{km}{h}$', ...
    'Interpreter', 'latex', 'Location','northeast')
%title('Dinamica del sistema controllato')

% con questi valori di k1 e k2 il bacino di attrazione è aumentato dato che
% le traiettorie convergono per tutte le condizioni inziali considerate,
% cosa che prima non accadeva

