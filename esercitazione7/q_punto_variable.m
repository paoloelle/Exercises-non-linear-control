function dqdt = q_punto_variable(t, q)

global k l1 l2 p p_punto

% p = @(t) [0; sin(t)];
% p_punto = @(t) [0; cos(t)];

dqdt = J_q(q)^-1*(k*(G_q(q) - p(t)) + p_punto(t));

end

