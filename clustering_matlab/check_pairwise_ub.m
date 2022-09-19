function [ml_violated, cl_violated] = check_pairwise_ub(X, ML, CL, tol)


    n_ml = size(ML, 1);
    n_cl = size(CL, 1);
    k = size(X, 2);
    
    ml_violated = 0;
    
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        fprintf('ML %d %d\n', i, j);
        if norm(X(i, :) - X(j, :)) > tol
            ml_violated = 1;
            break
        end
    end
    
    cl_violated = 0;

    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        fprintf('CL %d %d\n', i, j);
        if sum(abs(X(i, :) + X(j, :)) - 1 <= tol) < k
            cl_violated = 1;
            break
        end
    end

end