function [result] = call_solve_vector_lifting_root(data, card, params)

    rng(1727);
    disp(size(data))
    disp(card)
    disp(params)
    
    k = size(card, 1);
    result = solve_vector_lifting_child(data, card, [], [], cell(1, k), inf, [], params);
    result.best_B_cell = update_empty_B_cell(result.best_B_cell);
    if result.cp_flag == 1
        result.i_idx = -1;
        result.j_idx = -1;
    else
        [i_idx, j_idx, ~] = get_branching_pair(result.best_PI, card);
        result.i_idx = i_idx;
        result.j_idx = j_idx;
    end
    disp(result)
    
end
