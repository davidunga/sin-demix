function [comps, err] = stable_demix(v, Fs, w1, w2)
% A stable version of demix(). WIP.

arguments
    v
    Fs
    w1
    w2
end

v = v(:)';

% demix as usual:
[comps_full,opt] = linear_demix(v, Fs, [0, w1, w2]);

% susbtract DC, demix again:
comps_noDc = linear_demix(v - comps_full(1,:), Fs, [0, w1, w2]);

% final components: treat DC component of "no DC" demix as error, and
% therefore substract if from the DC component. Take sin components from
% "no DC" demix as-is:
comps = [comps_full(1,:) - comps_noDc(1,:); comps_noDc(2:3,:)];

% calc error:
vhat = sum(comps,1);
ii = opt.r:(length(v)-opt.r);
err = mean((vhat(ii)-v(ii)).^2/var(v(ii)));
