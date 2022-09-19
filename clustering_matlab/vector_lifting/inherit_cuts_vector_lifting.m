function [new_Bcell] = inherit_cuts_vector_lifting(ML_graph, init_B_cell, m, k)
    
    [bins_v , ~] = conncomp(ML_graph);
    
    %disp(bins_v)
    
    new_Bcell = cell(1, k);
    for j=1:k
        
        init_B_cell_j = init_B_cell{j};
        n_ineq_j = size(init_B_cell_j, 2);
        counter = 1;
        temp_Bcell = cell(1, n_ineq_j);
        for c=1:n_ineq_j
            
            [id_i, id_j, v] = find(init_B_cell_j{c});
            d = size(id_i, 1);
            
            for t=1:d
                
                id_i(t) = bins_v(id_i(t)-1);
                id_j(t) = bins_v(id_j(t)-1);
                
            end
        
            temp_Bcell{counter} = sparse(id_i+1, id_j+1, v, m+1, m+1);
            counter = counter + 1;
            
        end
        
        n_ineq_j = counter - 1;
        new_Bcell{j} = temp_Bcell(1:n_ineq_j);
        
    end
    
end