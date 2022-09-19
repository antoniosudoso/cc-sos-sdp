function out_B_cell = update_empty_B_cell(B_cell)

    k = size(B_cell, 2);
    out_B_cell = cell(1, k);
    for j=1:k
        if isempty(B_cell{j})
            out_B_cell{j} = {};
        else
            out_B_cell{j} = B_cell{j};
        end
    end

end