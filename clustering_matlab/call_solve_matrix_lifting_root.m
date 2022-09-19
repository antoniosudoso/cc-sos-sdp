function [result] = call_solve_matrix_lifting_root(data, card, params)

    rng(1727);
    disp(size(data))
    disp(card)
    disp(params)

    result = solve_matrix_lifting_child(data, card, [], [], cell(0), inf, [], params);
    if result.cp_flag == 1
        result.i_idx = -1;
        result.j_idx = -1;
    else
        [i_idx, j_idx, ~] = get_branching_pair(result.best_Z, card);
        result.i_idx = i_idx;
        result.j_idx = j_idx;
    end
    disp(result)
    
end
