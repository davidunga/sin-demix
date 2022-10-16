function plot_gt_and_prediction(t, comps, comps_hat)
% INPUTS:
% comps - components matrix
% t - time vec, length = size(comps,2)
% is_gt - components are of ground truth

comp_names = ["Signal", "DC", "Sin1", "Sin2"];

gt_colors = hsv(5);
gt_options = struct(LineWidth = 4, LineStyle = '-');

hat_colors = .2 * ones([5,3]);
hat_options = struct(LineWidth = 2, LineStyle = ':');

comps = [sum(comps,1);comps];
comps_hat = [sum(comps_hat,1);comps_hat];

figure();
tiledlayout(4,1);

for i = 1 : 4
    nexttile();
    hold on;
    plot(t, comps(i,:), 'color', gt_colors(i,:), gt_options);
    plot(t, comps_hat(i,:), 'color', hat_colors(i,:), hat_options);
    xlabel('t');
    title(comp_names(i));
end