function [x] = jacobi_2(G, x ,tol)   

    nG = length(x);
    b = zeros(nG,1);
    for i = 1:nG
        if x(i) ~= -1
            b(:,1) = b(:,1) - G(:,i)*x(i);          
        end
    end

    %Elimina Linhas e colunas
    G(x ~= -1, :) = [];
    G(:,x ~= -1) = [];
    b(x ~= -1, :) = [];
    n = length(x(x == -1));
    xk = 0.1*ones(n,1);
    xkp = xk;
    k = 1;
    while true
        for j = 1 : n
         xkp(j) = ((b(j) - G(j,[1:j-1,j+1:n]) * xk([1:j-1,j+1:n])) / G(j,j));
        end
        k = k + 1;
        if norm(xkp - xk, 2)/norm(xkp,2) < tol
            break
        end
        
        xk = xkp;
    end
    
    for i = 1:nG
        if x(i) == -1
            x(i) = xkp(1);
            xkp(1) = [];
        end
    end  
end