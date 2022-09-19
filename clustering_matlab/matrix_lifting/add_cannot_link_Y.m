function [B_cell_X_cl, l_X_cl, A_cell_Z_cl, b_X_cl] = add_cannot_link_Y(n, k, CL)
    
    n_cl = size(CL, 1);
    B_cell_X_cl = cell(1, k*n_cl);
    A_cell_Z_cl = cell(1, n_cl);
    l_X_cl = -ones(k*n_cl, 1);
    b_X_cl = zeros(n_cl, 1);
    counter = 1;
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        for h=1:k
            B_cell_X_cl{counter} = sparse([i+k, h, j+k, h], [h, i+k, h, j+k], [-0.5, -0.5, -0.5, -0.5], n+k, n+k);
            counter = counter + 1;
        end
        A_cell_Z_cl{c} = sparse([i+k, j+k], [j+k, i+k], [0.5, 0.5], n+k, n+k);
    end
end