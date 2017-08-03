clear
clc
load('petColor.mat');

qtdT = size(t,2);
qtdP = size(p,2);

[A, adj] = matEl(p, t, qtdP, qtdT);

b = cell(qtdT,1);
for i = 1:qtdT
    b{i} = ones(3,1);
end

[listaCor, RHS] = color(adj, t, qtdT, qtdP, b);
printColor(listaCor, p);