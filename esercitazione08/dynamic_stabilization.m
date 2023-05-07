function dxdt = dynamic_stabilization(t, x0)

global a k rho

x = x0(1); % condizione inziale stato
a_hat = x0(2); % condizione inziale stimatore

u = -a_hat*x - k*x; % legge di controllo

xd = a*x + u; % la a la metto per far evolvere il sistema anche se è incognita

a_hat_d = rho*x^2;

dxdt = [xd; a_hat_d];

end

