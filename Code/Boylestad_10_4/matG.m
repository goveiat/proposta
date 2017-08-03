function G = matG(C, t, qtdP, qtdT)
    G = sparse(qtdP,qtdP);
    for el=1:qtdT
        for i=1:3
            for j=1:3
                x = t(i,el);
                y = t(j,el);
                G(x,y) = G(x,y) + C{el}(i,j);
            end
        end
    end
end