clc
clear

% ebe | ebe_cor | ebe_cor | ebe_cor_par | ebe_cor_spmd | ebe_ot
alg = 'ebe_ot';
ref = 1;

load('pdetool/geometry.mat');
load('pdetool/coords.mat');

g = decsg(gd, sf, ns);

%Criação do modelo e da malha
m = createpde;
geometryFromEdges(m,g);
generateMesh(m, 'Hmax',ref);

%Obtenção dos pontos da malha
[p,e,t] = meshToPet(m.Mesh);

tic;

qtdT = size(t,2);
qtdP = size(p,2);

%Obtenção das matrizes elementares de coeficientes e alista de adj
[A, adj, map] = matEl(p, t, qtdP, qtdT);

%obtém as condições de x
cont = defContorno(p, qtdP);

%aplica as condições de contorno
[A, b] = setContorno_ebe(A, cont, adj,qtdT);


switch alg
    case 'ebe'
        [x, error, iter, flag, rhs] = cg_ebe(A, b, t, qtdP, qtdT, 1000, 0.001);
    case 'ebe_ot'
        [x, error, iter, flag, rhs] = cg_ebe_ot(A, b, t, qtdP, qtdT, 1000, 0.001, map);        
    case 'ebe_cor'
        startTimeCor = tic;
        listaCor = color(adj, t, qtdT);
        T = toc(startTimeCor);
        [x, error, iter, flag, tempo] = cg_ebe_cor(A, b,listaCor, qtdP, 1000, 0.001,m); 
    case 'ebe_cor_par'
        startTimeCor = tic;
        [A,MAP, RHS, qtdCor] = colorMap(adj, t, qtdT, qtdP, b, A);
        T = toc(startTimeCor);
        [x, error, iter, flag, tempo] = cg_ebe_par(A, MAP, RHS, qtdCor, qtdP, 1000, 0.001, m);        
    case 'ebe_cor_spmd'
        startTimeCor = tic;
        [A,MAP, RHS, qtdCor] = colorMap(adj, t, qtdT, qtdP, b, A);
        T = toc(startTimeCor);
        [x, error, iter, flag, tempo] = cg_ebe_spmd(A, MAP, RHS, qtdCor, qtdP, 1000, 0.001, m);               
end


toc;

% Impressões
COLOR =     0;
MESH =      0;
POT2D =     0;
POT3D =     1;
VECCAMPO =  0;
MATCAMPO =  0;


if COLOR
    figCor = printColor(listaCor, p);
    exportPDF(figCor, alg, ref, 'color')   
end



if MESH
    figMesh = figure;
    pdeplot(m)
    hold on
    scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')
    exportPDF(figMesh, alg, ref, 'mesh')     
end



if POT2D
    fig2DPot = figure;
    cx = pdetool_coords_x;
    cy = pdetool_coords_y;
    s = cell(8,1);
    for k = 1:8
        s{k} = p(:, p(1,:) <= cx{k}(2) & p(1,:) >= cx{k}(1) & p(2,:) >= cy{k}(1) & p(2,:) <= cy{k}(2));
    end
    S = [s{1} s{2} s{3} s{4} s{5} s{6} s{7} s{8}];
    pdeplot(m,'XYData',x','ColorMap','jet')
    hold on
    scatter(S(1,:),S(2,:), 15)
    exportPDF(fig2DPot, alg, ref, '2DPot')     
end



if VECCAMPO
    figVecCampo = figure;
    [ux,uy] = pdegrad(p,t,x'); 
    pdeplot(m,'FlowData',[ux;uy])
    hold on
    scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')
    exportPDF(figVecCampo, alg, ref, 'vecCampo')    
end



if MATCAMPO
    figMatizCampo = figure;
    pdeplot(m,'XYData',[ux;uy])
    hold on
    scatter(pts(1,:),pts(2,:), 5, 'FaceColor','r')
    exportPDF(figMatizCampo, alg, ref, 'matizCampo')
end



if POT3D
    fig3DPot = figure;
    pdeplot(m,'XYData',x,'ZData',x, 'ColorMap','jet')
    exportPDF(fig3DPot, alg, ref, '3DPot')
end









