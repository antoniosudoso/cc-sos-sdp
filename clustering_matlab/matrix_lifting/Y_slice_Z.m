function [A_rowsum, A_diag] = Y_slice_Z(n, k, C)

    A_rowsum = cell(1, n);
    A_diag = cell(1, n);
    C_vec = diag(C)';
    
    c = 1;
    for i = (k+1):(n+k)
        A_rowsum{c} = sparse([repelem(i, n), k+1:n+k], [k+1:n+k, repelem(i, n)], 0.5, n+k, n+k);
        %disp(full(A_rowsum{c}))
        coeff = [repelem(-0.5, k) ./ C_vec, repelem(-0.5, k) ./ C_vec, 1];
        A_diag{c} = sparse([repelem(i, k), 1:k, i], [1:k, repelem(i, k), i], coeff, n+k, n+k);
        %disp(full(A_diag{c}))
        c = c + 1;
    end
    
end