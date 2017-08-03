function r = setRHS(e,b)
    map = e([1 2 3]);
    rh =  b{e(5)};
    r = [map;rh];
end
