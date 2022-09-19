function [A_cell_pi_ml, b_pi_ml, A_cell_PI_ml, b_PI_ml] = add_must_link_P(n, ML)
    
    n_ml = size(ML, 1);
    A_cell_pi_ml = cell(1, n_ml);
    b_pi_ml = zeros(n_ml, 1);
    A_cell_PI_ml = cell(1, n*n_ml);
    b_PI_ml = zeros(n*n_ml, 1);
    counter_pi = 1;
    counter_PI = 1;
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        A_cell_pi_ml{counter_pi} = sparse([i+1, 1, j+1, 1], [1, i+1, 1, j+1], [0.5, 0.5, -0.5, -0.5], n+1, n+1);
        counter_pi = counter_pi + 1;
        for h=1:n
            A_cell_PI_ml{counter_PI} = sparse([i+1, h+1, j+1, h+1], [h+1, i+1, h+1, j+1], [0.5, 0.5, -0.5, -0.5], n+1, n+1);
            counter_PI = counter_PI + 1;
        end
    end
end