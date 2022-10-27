function evaluate_prediction(comps_gt, comps_hat, r)
% Evaluate prediction (comps_hat) relative to ground-truth (comps_gt)

%r = find(comps_hat(1,:)~=comps_hat(1,1),1,'first');
comps_hat = comps_hat(:,r:end-r);
comps_gt = comps_gt(:,r:end-r);

comp_names = ["DC", "Sin1", "Sin2"];

disp('Abs Relative Error:');
relerror(sum(comps_gt,1), sum(comps_hat,1), 'signal');
for i = 1 : length(comp_names)
    relerror(comps_gt(i, :), comps_hat(i, :), comp_names(i));
end

function err = relerror(gt, hat, name)
amp = std(gt,[],2)';

err = (hat - gt) ./ amp;
fprintf(' %-8s Mean=%2.2f Median=%2.3f Max=%2.3f\n', name, mean(abs(err)), median(abs(err)), max(abs(err)));