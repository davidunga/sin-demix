function plot_params(params, t, is_gt, lgnd)

colors = hsv(5);

if is_gt
    options.LineWidth = 4;
    options.LineStyle = '-';
else
    options.LineWidth = 2;
    options.LineStyle = ':';
    colors = .2 * ones(size(colors));
end

hold on;
paramNames = sort(setdiff(fieldnames(params), {'w1', 'w2'}));
%paramNames = {'p1','p2'};
for i = 1 : length(paramNames)
    xx = params.(paramNames{i});
    %xx = (xx - mean(xx)) / std(xx);
    plot(t, xx, 'color', colors(i,:), options);
end

if exist('lgnd', 'var') && lgnd
    legend(paramNames);
end