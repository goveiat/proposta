clc

%Definição da Geometria do Espaço e do Capacitor
% geo = [3;4;xcoords;ycoords]
espaco = [3 4 0 16.0794 16.0794 0 17.08 17.08 0 0]';
capacitor = [3 4 8 8.0794 8.0794 8 11.08 11.08 6 6]';
g = decsg([espaco capacitor]);

model = createpde;
geometryFromEdges(model,g);
generateMesh(model);


%Criação da malha
[p,ed,t] = initmesh(g);

%Refinamento
for i = 1:2
    [p,ed,t] = refinemesh(g,p,ed,t);
end

%Relação Refinamento nº de Pontos
%   1   2.433       pontos
%   2   9.649       pontos
%   3   38.433      pontos
%   4   153.409     pontos
%   5   612.993     pontos
%   6   2.450.689   pontos

%Obtenção dos pontos das placas
placaD = p(:, p(1,:) == 8.0794 & p(2,:) >= 6 & p(2,:) <= 11.08);
placaE = p(:, abs(p(1,:)-8) <= eps(8)  & p(2,:) >= 6 & p(2,:) <= 11.08);
interior = p(:, p(1,:) > 8 & p(1,:) < 8.0794 & p(2,:) >= 6 & p(2,:) <= 11.08);

%Impressão da malha
figure
pdeplot(model)
hold on
scatter([placaE(1,:) interior(1,:) placaD(1,:)],[placaE(2,:) interior(2,:) placaD(2,:)], 3, 'MarkerEdgeColor','k', 'MarkerFaceColor','k')


qtdT = size(t,2);
qtdP = size(p,2);

C = zeros(qtdT, 3, 3);
for e=1:qtdT  
    coords = p(:,t([1 2 3],e));
    area = 0.5 * det([coords' [1;1;1]]);
    if area < 0
        coords = flipud(coords);
        area = 0.5 * det([coords' [1;1;1]]);
    end

    P = coords(2, [2 3 1]) - coords(2, [3 1 2]);
    Q = coords(1, [2 3 1]) - coords(1, [3 1 2]);
    for i = 1:3
        for j= 1:3
            C(e,i,j) = (P(i)*P(j) + Q(i)*Q(j))/(4*area);
        end
    end
end

G = zeros(qtdP);
for e=1:qtdT
    for i=1:3
        for j=1:3
            x = t(i,e);
            y = t(j,e);
            G(x,y) = G(x,y) + C(e,i,j);
        end
    end
end

vEspaco = -ones(1, qtdP);
vPlacaE = 49*ismember(p,placaE);
vPlacaD = ismember(p,placaD);

% Array com condições de contorno
contorno = vEspaco + vPlacaE(1,:) + vPlacaD(1,:);



F = find(contorno == -1);
P = find(contorno ~= -1);
Vp = contorno(P)';
Gff = G(F, F);
Gfp = G(F, P);
b = -Gfp*Vp;
Vf = Gff\b;
contorno(contorno == -1) = Vf;

x = p(1,:);
y = p(2,:);




figure
pdesurf(p,t,contorno');colorbar;



