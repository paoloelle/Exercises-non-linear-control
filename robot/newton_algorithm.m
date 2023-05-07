function dqdt = newton_algorithm(t, q)

global k p pd

dqdt = J_q(q)^-1*(k*(G_q(q) - p(t)) + pd(t));

