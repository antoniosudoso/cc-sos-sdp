function A = Y_zeros_cell(m, n)

    % m constraints

    A = cell(1, m);
 
    c = 1;
    for i = 1:m
        A{c} = sparse(n, n);
        c = c + 1;
    end
    
end