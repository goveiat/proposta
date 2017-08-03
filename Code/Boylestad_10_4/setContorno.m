function [Gff, b] = setContorno(G, x)
    F = find(x == -1);
    P = find(x ~= -1);
    Vp = x(P)';
    Gff = G(F, F);
    Gfp = G(F, P);
    b = -Gfp*Vp;
end