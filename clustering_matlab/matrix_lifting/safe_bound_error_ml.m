function [LB0, LB] = safe_bound_error_ml(blk, At, C, b, y, Z2, Bt, y2, l, Y, max_card)
  
  % Safe bounds through error bounds
  
  mu = 1.1;
  Aty = sdpnalAtyfun(blk, At, y);
  Znew = ops(C, '-', Aty);
  if ~isempty(Bt)
    Bty = sdpnalAtyfun(blk, Bt, y2);
    Znew = ops(Znew, '-', Bty);
  end
  if ~isempty(Z2)
     Znew = ops(Znew, '-', Z2); 
  end

  eigtmp = eig(full(Znew{1}));
  idx = find(eigtmp < -1e-7); 
  Xbar = min(mu * max(eig(full(Y{1}))), max_card + 1);
  LB0 = b'*y + l'*y2;
  pert = Xbar * sum(eigtmp(idx));
  LB = LB0 + pert; 
  %fprintf('PP safe lower bound = %10.9e \n', LB);
  
end