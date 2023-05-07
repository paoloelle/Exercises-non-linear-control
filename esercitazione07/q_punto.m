function dqdt = q_punto(t, q)
    
    global k p p_punto

    dqdt = J_q(q)^-1*(k*(G_q(q) - p) + p_punto);

end


