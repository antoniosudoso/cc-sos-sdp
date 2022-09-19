function [new_Bcell] = inherit_cuts_matrix_lifting(ML_graph, init_B_cell, m, k)
    
    n_ineq = size(init_B_cell, 2);
    
    new_Bcell = cell(1, n_ineq);
    [bins_v , ~] = conncomp(ML_graph);
    
    %disp(bins_v)
    
    counter = 1;
    for c=1:n_ineq
        
        [id_i, id_j, v] = find(init_B_cell{c});
        
        d = size(id_i, 1);
        
        for t=1:d
            
            id_i(t) = bins_v(id_i(t)-k);
            id_j(t) = bins_v(id_j(t)-k);

        end
        
        new_Bcell{counter} = sparse(id_i+k, id_j+k, v, m+k, m+k);
        counter = counter + 1;
                
    end
    
    n_ineq = counter - 1;
    new_Bcell = new_Bcell(1:n_ineq);
    %disp(new_Bcell);
    
end