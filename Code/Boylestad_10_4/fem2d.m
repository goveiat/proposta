clc
clear

load('pdetool/geometry.mat');
load('pdetool/coords.mat');

g = decsg(gd, sf, ns);

%Criação do modelo e da malha
m = createpde;
geometryFromEdges(m,g);
generateMesh(m, 'Hmax',0.3);

%Obtenção dos pontos da malha
[p,e,t] = meshToPet(m.Mesh);

tic;

qtdT = size(t,2);
qtdP = size(p,2);

%Obtenção das matrizes elementares de coeficientes
C = matEl(p, t, qtdP, qtdT);

%Obtenção da matriz global de coeficientes
G = matG(C, t, qtdP, qtdT);

%obtém as ondições de x
[x, pts] = defContorno(p, qtdP);

%Aplica o x na matriz
[G, b] = setContorno(G, x);

%Direto
sol = G\b;


%CG
% chute = zeros(length(b),1);
% M = length(b);
% sol = cg(G, chute, b, eye(M), 1000, 0.001);



x(x == -1) = sol;

toc;

% %% Impressões
% %Malha
%  figure
%  pdeplot(m)
% hold on
% scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')
% 
%Distribuição de Potencial 
cx = pdetool_coords_x;
cy = pdetool_coords_y;
s = cell(8,1);
for k = 1:8
    s{k} = p(:, p(1,:) <= cx{k}(2) & p(1,:) >= cx{k}(1) & p(2,:) >= cy{k}(1) & p(2,:) <= cy{k}(2));
end
S = [s{1} s{2} s{3} s{4} s{5} s{6} s{7} s{8}];
fig = figure
pdeplot(m,'XYData',x','ColorMap','jet')
% hold on
% pdeplot(m)
hold on
scatter(S(1,:),S(2,:), 15, 'MarkerEdgeColor', 'k')
% 
% %Vetor Campo Elétrico
% figure
% [ux,uy] = pdegrad(p,t,x'); 
% pdeplot(m,'FlowData',[ux;uy])
% hold on
% scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')
% 
% figure
% pdeplot(m,'XYData',[ux;uy])
% hold on
% scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')

% figure
% pdeplot(m,'XYData',x,'ZData',x, 'ColorMap','jet')

set(fig,'Units','Inches');
    pos = get(fig,'Position');
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    nmFig = ['_0.2' '.pdf'];
    print(fig,nmFig,'-dpdf','-r0')








