function rhs = setRHS2(e,b,qtdN)
    rhs = sparse(qtdN,1); 
    map = e([1 2 3]);
    rhs(map) =  b{e(5)};
end