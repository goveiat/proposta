clc

%Definição da Geometria do Espaço e do Capacitor
espaco = [3 4 0 16.0794 16.0794 0 17.08 17.08 0 0]';
%espaco = [3 4 6 10 10 6 13 13 4 4]';
capacitor = [3 4 8 8.0794 8.0794 8 11.08 11.08 6 6]';
g = decsg([espaco, capacitor]);

%Criação do modelo e da malha
m = createpde;
geometryFromEdges(m,g);
generateMesh(m,'Hmax',0.6);

%Obtenção dos pontos da malha
[p,ed,t] = meshToPet(m.Mesh);


qtdT = size(t,2);
qtdP = size(p,2);

%Obtenção das matrizes elementares de coeficientes
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


%Obtenção da matriz global de coeficientes
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

%Obtenção dos pontos das placas
placaD = p(:, p(1,:) == 8.0794 & p(2,:) >= 6 & p(2,:) <= 11.08);
placaE = p(:, abs(p(1,:)-8) <= eps(8)  & p(2,:) >= 6 & p(2,:) <= 11.08);
interior = p(:, p(1,:) > 8 & p(1,:) < 8.0794 & p(2,:) >= 6 & p(2,:) <= 11.08);

vPlacaE = 49*ismember(p,placaE);
vPlacaD = ismember(p,placaD);
vEspaco = -ones(1, qtdP);
contorno = vEspaco + vPlacaE(1,:) + vPlacaD(1,:);

pts = [placaE placaD];

%Matriz de banda
%contorno = mBanda(G, contorno);

%Jacobi
contorno= jacobi(G, contorno, 0.0001);


x = p(1,:);
y = p(2,:);

%% Impressões
%Malha
figure
pdeplot(m)
hold on
scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')

%Distribuição de Potencial 
figure
pdeplot(m,'XYData',contorno')
hold on
scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')

%Vetor Campo Elétrico
figure
[ux,uy] = pdegrad(p,t,-contorno'); 
pdeplot(m,'FlowData',[ux;uy])
hold on
scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')

figure
pdeplot(m,'XYData',[ux;uy])
hold on
scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')









