function [A_cell_X_ml, b_X_ml, A_cell_Z_ml, b_Z_ml] = add_must_link_Y(n, k, ML)
    
    n_ml = size(ML, 1);
    A_cell_X_ml = cell(1, k*n_ml);
    b_X_ml = zeros(k*n_ml, 1);
    A_cell_Z_ml = cell(1, n*n_ml);
    b_Z_ml = zeros(n*n_ml, 1);
    counter_X = 1;
    counter_Z = 1;
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        for h=1:k
            A_cell_X_ml{counter_X} = sparse([i+k, h, j+k, h], [h, i+k, h, j+k], [0.5, 0.5, -0.5, -0.5], n+k, n+k);
            counter_X = counter_X + 1;
        end
        for h=1:n
            A_cell_Z_ml{counter_Z} = sparse([i+k, h+k, j+k, h+k], [h+k, i+k, h+k, j+k], [0.5, 0.5, -0.5, -0.5], n+k, n+k);
            counter_Z = counter_Z + 1;
        end
    end
end