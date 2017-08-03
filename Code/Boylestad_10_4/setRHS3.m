function rhs = setRHS3(t,b,qtdN)
    rhs = zeros(qtdN,1); 
    for e = t
        map = e([1 2 3]);
        rhs(map) =  b{e(5)};        
    end
end