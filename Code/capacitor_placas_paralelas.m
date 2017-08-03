clear
clc
close all

r_tam = [2 2];
r_div = [20 20];
r_passo = r_tam ./ r_div;
n_areaElf = prod(r_passo)/2;
n_qtdElf = prod(r_div)*2;
r_qtdNos= [r_div(1)+1 r_div(2)+1];
n_qtdNos= prod(r_qtdNos);
m_nos = zeros(n_qtdNos,2);
m_elf = zeros(n_qtdElf,3);
r_v = -ones(1, n_qtdNos);
m_v = reshape(r_v, r_qtdNos);
r_b = zeros(1, n_qtdNos);

k = 1;
for i = 1:r_qtdNos(1)
    for j = 1:r_qtdNos(2)
        m_nos(k, :) = [(j-1)*r_passo(1) (i-1)*r_passo(2)];
        k=k+1;
    end
end


for e = 1:n_qtdElf
    k = ceil(e/(2*r_div(1)))-1;
    if rem(e,2)==1
        m_elf(e,1)=(e+1)/2+r_div(1)+1+k;
        m_elf(e,2)=(e+1)/2+k;
        m_elf(e,3)=(e+1)/2+1+k;
    else
        m_elf(e,1)=e/2+r_div(1)+1+k;
        m_elf(e,2)=e/2+1+k;
        m_elf(e,3)=e/2+r_div(1)+2+k;
    end
end


figure
subplot(2,2,1);
TR = triangulation(m_elf,m_nos);
hold on
triplot(TR, 'Color','k')
rectangle('Position',[0.9 0.5 0.1 1], 'FaceColor','r')
rectangle('Position',[1.1 0.5 0.1 1], 'FaceColor', 'b')
ptsP1 = m_nos(m_nos(:,1) >= 0.9 & m_nos(:,1) <= 1 & m_nos(:,2) >= 0.5 & m_nos(:,2) <= 1.5, :);
maxP1 = max(ptsP1,[],1);
ptsP2 = m_nos(m_nos(:,1) >= 1.1 & m_nos(:,1) <= 1.2 & m_nos(:,2) >= 0.5 & m_nos(:,2) <= 1.5, :);
minP2 = min(ptsP2,[],1);
if isempty(minP2) || isempty(maxP1)
    error('Triangulaçao mal formada');
end
linhaP1 = ptsP1(ptsP1(:,1) == maxP1(1),:);
linhaP2 = ptsP2(ptsP2(:,1) == minP2(1),:);
pts = [linhaP1; linhaP2];
scatter(pts(:,1),pts(:,2),15, 'FaceColor', 'k')
xlabel('x');
ylabel('y');
title ('Triangulação do Espaço e Pontos');


r_v(find(ismember(m_nos,linhaP1,'rows'))) = 10;
r_v(find(ismember(m_nos,linhaP2,'rows'))) = 0;


C = zeros(n_qtdElf, 3, 3);
for e=1:n_qtdElf
    P = m_nos(m_elf(e,[2 3 1]),2)-m_nos(m_elf(e,[3 1 2]),2);
    Q = m_nos(m_elf(e,[2 3 1]),1)-m_nos(m_elf(e,[3 1 2]),1);
%     n_areaElf = 0.5 * (P(2)*Q(3) - P(3)*Q(2));
    for i = 1:3
        for j= 1:3
            C(e,i,j) = (P(i)*P(j) + Q(i)*Q(j))/(4*n_areaElf);
        end
    end
end

G = zeros(n_qtdNos);
for e=1:n_qtdElf
    for i=1:3
        for j=1:3
            n_x = m_elf(e,i);
            n_y = m_elf(e,j);
            G(n_x,n_y) = G(n_x,n_y) + C(e,i,j);
        end
    end
end

f = find(r_v == -1);
p = find(r_v ~= -1);
Vp = r_v(p)';
Gff = G(f, f);
Gfp = G(f, p);
b = -Gfp*Vp;
Vf = Gff\b;
r_v(r_v == -1) = Vf;
m_v = reshape(r_v, r_qtdNos)';

x = 0:r_passo(1):r_tam(1);
y = 0:r_passo(2):r_tam(2);

%figure
subplot(2,2,2);
pcolor(x,y,m_v);shading interp; colorbar;
hold on
scatter(pts(:,1),pts(:,2), 3, 'FaceColor','r')
xlabel('x');
ylabel('y');
title ('Distribuição do Potencial (V)');

%figure
subplot(2,2,3);
[Ex,Ey]=gradient(m_v);
Ex=-Ex/r_passo(1);
Ey=-Ey/r_passo(2);

xlabel('x');
ylabel('y');
[C h] = contour(x,y,m_v);
clabel(C,h);
hold on
quiver(x,y,Ex,Ey)
scatter(pts(:,1),pts(:,2), 5, 'FaceColor','r')
title(' Distribuição do Campo Elétrico e Linhas de Potencial');

%figure
subplot(2,2,4);
surf(y, x,Ey, Ex);shading interp; colorbar;
hold on
scatter(pts(:,1),pts(:,2), 3, 'FaceColor','r')
xlabel(' x');
ylabel(' y');
title(' Intensidade do Campo Elétrico (V/m)');
view([0 90]);

figure
streamslice(x, y, Ex, Ey);
xlabel(' x');
ylabel(' y');
title(' Linhas do Campo Elétrico');
hold on
scatter(pts(:,1),pts(:,2), 3, 'FaceColor','r')



