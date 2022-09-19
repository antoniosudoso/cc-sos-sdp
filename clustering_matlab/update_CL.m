function [new_CL] = update_CL(ML_graph, CL)

    n_cl = size(CL, 1);
    new_CL = zeros(n_cl, 2);
    [bins_v , ~] = conncomp(ML_graph);
    for c=1:n_cl
        id_comp_i = bins_v(CL(c, 1));
        id_comp_j = bins_v(CL(c, 2));
        if id_comp_i < id_comp_j
            new_CL(c, :) = [id_comp_i, id_comp_j];
        else
            new_CL(c, :) = [id_comp_j, id_comp_i];
        end
    end
    new_CL = unique(new_CL, 'rows');
    
end