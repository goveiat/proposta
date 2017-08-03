clc
clear

%Definição da Geometria do Espaço e do Capacitor
espaco = [3 4 0 16.0794 16.0794 0 17.08 17.08 0 0]';
%espaco = [3 4 6 10 10 6 13 13 4 4]';
capacitor = [3 4 8 8.0794 8.0794 8 11.08 11.08 6 6]';
g = decsg([espaco, capacitor]);

%Criação do modelo e da malha
m = createpde;
geometryFromEdges(m,g);
generateMesh(m,'Hmax',0.3);

% applyBoundaryCondition(m,'Edge',8,'u',10);
% applyBoundaryCondition(m,'Edge',3,'u',0);

% specifyCoefficients(m,'m',0,'d',0,'c',1,'a',0,'f',0);
% generateMesh(m,'Hmax',0.3);
% results = solvepde(m);
% u = results.NodalSolution;
% pdeplot(model,'XYData',u,'ZData',u)
% view(-23,8)

figure
pdeplot(m);
% figure
% pdegplot(m,'EdgeLabels','on');