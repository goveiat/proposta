function [fig] = printColor(listaCor, p)
    n = length(listaCor);

    fig = figure;
    rng('default')
    for i = 1:n
        x = [];
        y = [];
        for ee = listaCor{i}
            no = p(:, ee([1 2 3], 1)')';
            x = [x no(:,1)];
            y = [y no(:,2)];
        end
    patch(x,y,rand(1,3))
    rng(i)
    end
end