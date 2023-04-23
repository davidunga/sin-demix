function [model, inlier_ixs] = ransac_line(pts, inlier_margin)
rng(0);

fit_fcn = @(pts) polyfit(pts(:,1),pts(:,2),1);
eval_fcn = @(model, pts) sum((pts(:, 2) - polyval(model, pts(:,1))).^2,2);

[~, inlier_ixs] = ransac(pts, fit_fcn, eval_fcn, 2, inlier_margin^2);
model = polyfit(pts(inlier_ixs,1), pts(inlier_ixs,2),1);