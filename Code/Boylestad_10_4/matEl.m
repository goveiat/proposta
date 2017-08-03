function [C, adj, map] = matEl(p, t, qtdN, qtdT)

    C = cell(qtdT, 1);
    adj = cell(qtdN, 1);
    for e=1:qtdT  
        nT = t([1 2 3],e);
        map(:,e) = nT;
        k = 1;
        for n = nT'
            adj{n}([1 2],end+1) = [e; k];
            k = k+1;
        end
        
    
        coords = p(:,nT);
        area = 0.5 * det([coords' [1;1;1]]);
        if area < 0
            coords = flipud(coords);
            area = 0.5 * det([coords' [1;1;1]]);
        end

        P = coords(2, [2 3 1]) - coords(2, [3 1 2]);
        Q = coords(1, [2 3 1]) - coords(1, [3 1 2]);
        temp = zeros(3);
        for i = 1:3
            for j= 1:3
                temp(i,j) = (P(i)*P(j) + Q(i)*Q(j))/(4*area);
            end
        end
        C{e} = temp;
    end
end