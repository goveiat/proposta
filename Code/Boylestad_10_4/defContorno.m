function [contorno, pts] = defContorno(p, qtdP)
    placaE = p(:, (p(1,:) == 7.9603 | p(1,:) == 7.9206) & (p(2,:) == 5.46 | p(2,:) == 10.54));
    placaD = p(:, (p(1,:) == 8.0397 | p(1,:) == 8.0794) & (p(2,:) == 5.46 | p(2,:) == 10.54));

    vPlacaE = 49*ismember(p,placaE);
    vPlacaD = ismember(p,placaD);
    vEspaco = -ones(1, qtdP);
    contorno = vEspaco + vPlacaE(1,:) + vPlacaD(1,:);

    pts = [placaD placaE];
end