function [comps, err] = stable_demix(v, Fs, w1, w2)

arguments
    v
    Fs
    w1
    w2
end

v = v(:)';

[comps_full,opt] = linear_demix(v, Fs, [0, w1, w2]);

comps_noDc = linear_demix(v - comps_full(1,:), Fs, [0, w1, w2]);
comps = [comps_full(1,:) - comps_noDc(1,:); comps_noDc(2:3,:)];

vhat = sum(comps,1);
ii = opt.r:(length(v)-opt.r);
err = (vhat(ii)-v(ii)).^2/var(v(ii));