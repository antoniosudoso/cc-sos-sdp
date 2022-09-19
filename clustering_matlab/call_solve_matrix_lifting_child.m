function [result] = call_solve_matrix_lifting_child(data, card, ML, CL, B_cell, global_ub, global_X, params)

    rng(1727);

    disp('ML')
    disp(ML)
    disp('CL')
    disp(CL)

    result = solve_matrix_lifting_shrinking(data, card, ML, CL, B_cell, global_ub, global_X, params);
    result.best_ub = min(result.ub_list);
    if result.cp_flag == 1
        result.i_idx = -1;
        result.j_idx = -1;
    else
        [i_idx, j_idx, ~] = get_branching_pair(result.best_Z, card);
        result.i_idx = i_idx;
        result.j_idx = j_idx;
    end
    %disp(result)
    
end