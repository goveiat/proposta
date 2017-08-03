clc
%Definição da Geometria do Espaço e do Capacitor
E = [3 4 0 16.0794 16.0794 0 17.08 17.08 0 0]';
C = [3 4 8 8.0794 8.0794 8 11.08 11.08 6 6]';
gd = [E,C];
sf = 'E+C';
ns = char('E','C');
ns = ns';
g = decsg(gd,sf, ns);

%tool
%pderect([0 16 0 16])
%pderect([7.9206 7.9603 5.46 10.54])
%pderect([8.0397 8.0794 5.46 10.54])

%Criação do modelo e da malha
model = createpde;
geometryFromEdges(model,g);

figure
pdegplot(model,'EdgeLabels','on');
generateMesh(model,'Hmax',0.5);
[p,e,t] = meshToPet(model.Mesh);

applyBoundaryCondition(model,'Edge',8,'u',0);
applyBoundaryCondition(model,'Edge',3,'u',48);
u = assempde(model,1,0,0);

figure
pdeplot(model,'XYData',u)

figure
[ux,uy] = pdegrad(p,t,-u); 
pdeplot(model,'FlowData',[ux;uy])
axis equal

figure
pdeplot(model)
axis equal

figure
pdeplot(model,'XYData',[ux;uy], 'ColorMap','jet')
axis equal
