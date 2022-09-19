function [A_rowsum, A_diag] = Y_slice_Z_shrinking(m, k, C, e)

    % e = T^s 1_n (new size is m)

    A_rowsum = cell(1, m);
    A_diag = cell(1, m);
    C_vec = diag(C)';
    
    c = 1;
    for i = (k+1):(m+k)
        %A_rowsum{c} = sparse([repelem(i, m); k+1:m+k], [k+1:m+k; repelem(i, m)], 0.5 .* [e; e], m+k, m+k);
        A_rowsum{c} = sparse([repmat(i, m, 1); (k+1:m+k)'], [(k+1:m+k)'; repmat(i, m, 1)], [0.5 .* ones(m, 1) .* e; 0.5 .* ones(m, 1) .* e], m+k, m+k);
        %disp(full(A_rowsum{c}))
        coeff = [repelem(-0.5, k) ./ C_vec, repelem(-0.5, k) ./ C_vec, 1];
        A_diag{c} = sparse([repelem(i, k), 1:k, i], [1:k, repelem(i, k), i], coeff, m+k, m+k);
        %disp(full(A_diag{c}))
        c = c + 1;
    end
    
end