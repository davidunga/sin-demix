function plot_sines(comps, t, is_gt, is_sp, comp_names)
% INPUTS:
% comps - components matrix
% t - time vec, length = size(comps,2)
% is_gt - components are of ground truth
% is_sp - using subplot
% comp_names - component names

if ~exist('is_sp', 'var'), is_sp=false; end
if ~exist('comp_names', 'var'), comp_names=[]; end

colors = hsv(5);

if is_gt
    options.LineWidth = 4;
    options.LineStyle = '-';
else
    options.LineWidth = 2;
    options.LineStyle = ':';
    colors = .2 * ones(size(colors));
end


if is_sp
    nexttile(1);
    xlabel('t');
    if ~isempty(comp_names), title('signal'); end
end

hold on; plot(t, sum(comps,1), 'color', colors(1,:), options);
for i = 1 : 3
    if is_sp
        nexttile(i + 1);
        xlabel('t');
        if ~isempty(comp_names), title(comp_names{i}); end
    end
    hold on; plot(t, comps(i, :), 'color', colors(i + 1,:), options);
end

if ~is_sp && ~isempty(comp_names)
    legend(comp_names);
end
