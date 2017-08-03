function [w] = solvePar2(t,e,A,r)
  map = t([1 2 3],e)';
  M = diag(A{e});
  m = power(M,-1);    
  h = r(map);
  w = m.*h;
end