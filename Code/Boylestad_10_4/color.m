function [listaCor] = color(adj, t, qtdT)
t(4,:) = 0;
cores = [];
for el = 1:qtdT
    c = t([1 2 3],el)';
    viz = [];
    for cc = c
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
        listaCor{nova} = t(:,el);
    else
        min = 1000;
        for d = disponiveis
            s = size(listaCor{d},2);
            if s < min
                sel = d;
            end
        end
        %sel = min(disponiveis);
        pos = size(listaCor{sel},2)+1;
        t([4;5;6],el) = [sel el pos];
        listaCor{sel}(:,pos) = t(:,el);
    end
end


