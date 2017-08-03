function Q = dotPar(Ld, Pd)
    s = size(Ld,1);
    Q = zeros(1,s);
    for i=1:3:s
        j = 1+(i-1)/3;
        k = [i i+1 i+2];
        Q(:,k) = Ld(k,:)*Pd(:,j);
    end
end