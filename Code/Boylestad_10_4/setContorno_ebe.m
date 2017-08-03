function [A, b] = setContorno_ebe(A, cont, adj, qtdT)

b = cell(qtdT,1);
for e = 1:qtdT
    b{e} = zeros(3,1);
end

nos = find(cont ~= -1);
for no = nos
    val = cont(no);
    els = adj{no};
    qtd = size(els,2);
    for i = els
        e = i(1);
        ni = i(2);
        b{e} = b{e} - A{e}(:,ni)*val;
        b{e}(ni) = val/qtd;          
        A{e}(:,ni) = 0;
        A{e}(ni,:) = 0;
        A{e}(ni, ni) = 1/qtd;
    end
end