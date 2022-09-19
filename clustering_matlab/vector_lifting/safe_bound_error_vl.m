function [LB0, LB] = safe_bound_error_vl(blk, At, C, b, y, Z2, Bt, y2, l, Y, card)

  % Safe bounds through error bounds

  mu = 1.1;
  k = size(card, 1);
  sum_eig_list = zeros(k, 1);
  Xbar_list = zeros(k, 1);
  for j=1:k  
      Aty = sdpnalAtyfun(blk(j, :), At{j}, y);
      Znew = ops(C{j}, '-', Aty);
      if ~isempty(Bt)
        Bty = sdpnalAtyfun(blk(j, :), Bt{j}, y2);
        Znew = ops(Znew, '-', Bty);
      end
      if ~isempty(Z2{j})
         Znew = ops(Znew, '-', Z2{j}); 
      end
      Yj = Y{j};
      eigtmp = eig(full(Znew)); 
      idx = find(eigtmp < -1e-7); 
      sum_eig_list(j) = sum(eigtmp(idx));  
      Xbar_list(j) = min(mu * max(eig(full(Yj))), card(j,j) + 1);
  end
  
  LB0 = b'*y + l'*y2;
  pert = sum(sum_eig_list .* Xbar_list);
  LB = LB0 + pert; 
  %fprintf('PP safe lower bound = %10.9e \n', LB);
  
end