function [G,MAP, RHS, qtdCor] = colorMap(adj, t, qtdT, qtdP,b, A)
t(4,:) = 0;
cores = [];
RHS = zeros(qtdP,1);
for el = 1:qtdT
    map = t([1 2 3],el)';
    RHS(map) = RHS(map) + b{el};
    viz = [];
    for cc = map
        ad = adj{cc}(1,:);
        viz = [viz ad(ad ~= el)];
    end
    viz = unique(viz);
    corViz = t(4,viz);
    disponiveis = setdiff(cores, corViz);
    if isempty(disponiveis)
        nova = length(cores)+1;
        cores(1,end+1) = nova;
        t([4;5;6],el) = [nova el 1];
        qtdCor(nova) = 1;
        MAP{nova} = map';
        G{nova,1} = A{el};
    else
        min = 1000;
        for d = disponiveis
            s = qtdCor(d);
            if s < min
                sel = d;
            end
        end
        %sel = min(disponiveis);
        pos = qtdCor(sel)+1;
        t([4;5;6],el) = [sel el pos];
        qtdCor(sel) = qtdCor(sel)+1;        
        MAP{sel} = [MAP{sel}; map'];
        G{sel,pos} = A{el};
    end
end


