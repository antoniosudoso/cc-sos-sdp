function [B_cell_pi_cl, l_pi_cl, A_cell_PI_cl, b_PI_cl] = add_cannot_link_P(n, CL)
    
    n_cl = size(CL, 1);
    B_cell_pi_cl = cell(1, n_cl);
    A_cell_PI_cl = cell(1, n_cl);
    l_pi_cl = -ones(n_cl, 1);
    b_PI_cl = zeros(n_cl, 1);
    counter_pi = 1;
    counter_PI = 1;
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        B_cell_pi_cl{counter_pi} = sparse([i+1, 1, j+1, 1], [1, i+1, 1, j+1], [-0.5, -0.5, -0.5, -0.5], n+1, n+1);
        counter_pi = counter_pi + 1;
        A_cell_PI_cl{counter_PI} = sparse([i+1, j+1], [j+1, i+1], [0.5, 0.5], n+1, n+1);
        counter_PI = counter_PI + 1;
    end
end