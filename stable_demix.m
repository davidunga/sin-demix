function [comps, err, opt] = stable_demix(v, Fs, w1, w2, options)
% A stable version of demix(). WIP.

arguments
    v
    Fs
    w1
    w2
    options.mode = "naive"
end

disp(options.mode);

switch options.mode
    case "linear"
        demix_fnc = @linear_demix;
    case "wt"
        demix_fnc = @wt_demix;
    case "wt2"
        demix_fnc = @wt_demix2;
    case "wt3"
        demix_fnc = @wt_demix3;
    case "naive"
        demix_fnc = @naive_demix;
    otherwise
        error("Unknown demix mode %s", options.mode);
end

v = v(:)';

% demix as usual:
[comps_full,opt] = demix_fnc(v, Fs, [0, w1, w2]);

% susbtract DC, demix again:
comps_noDc = demix_fnc(v - comps_full(1,:), Fs, [0, w1, w2]);

% final components: treat DC component of "no DC" demix as error, and
% therefore substract if from the DC component. Take sin components from
% "no DC" demix as-is:
comps = [comps_full(1,:) - comps_noDc(1,:); comps_noDc(2:3,:)];

% calc error:
vhat = sum(comps,1);
ii = opt.r:(length(v)-opt.r);
err = mean((vhat(ii)-v(ii)).^2/var(v(ii)));
