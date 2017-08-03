function Z = solvePar(Ld, Rd)
    s = size(Ld,1);
    Z = zeros(1,s);
    for i=1:3:s
        j = 1+(i-1)/3;
        k = [i i+1 i+2];
        m = power(diag(Ld(k,:)),-1); 
        Z(:,k) = m.*Rd(:,j);
    end
end