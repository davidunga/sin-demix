function evaluate_prediction(comps_gt, comps_hat)
% Evaluate prediction (comps_hat) relative to ground-truth (comps_gt)

comp_names = ["DC", "Sin1", "Sin2"];

disp('Abs Relative Error:');
relerror(sum(comps_gt,1), sum(comps_hat,1), 'signal');
for i = 1 : length(comp_names)
    relerror(comps_gt(i, :), comps_hat(i, :), comp_names(i));
end

function err = relerror(gt, hat, name)
if strcmpi(name, 'dc')
    amp = mean(gt);
else
    amp = diff(prctile(gt, [1, 99]));
end

err = (hat - gt) ./ amp;
fprintf(' %-8s Median=%2.3f Max=%2.3f\n', name, median(abs(err)), max(abs(err)));