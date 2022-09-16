function plot_sines(comps, t, is_gt)

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
plot(t, sum(comps,1), 'color', colors(1,:), options);
for i = 1 : 3
    plot(t, comps(i, :), 'color', colors(i + 1,:), options);
end