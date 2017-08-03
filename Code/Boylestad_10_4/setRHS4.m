function rhs = setRHS4(b,qtdN,t)
      rhs = zeros(qtdN,1);
      s = size(b,2);
      for e = 1:s
          map = t([1 2 3],e)';
          rhs(map,1) = rhs(map,1) +  b(:,e);
      end
end