function [ws_refined,ph,margin,rest_ixs] = refine_sine_params(v, Fs, ws)

arguments
    v
    Fs
    ws
end

% ---
dA_THRESH = .00005;
MIN_REST_DUR = .05;
W_REFINE_TOL = .01;
COI_THRESH = 50;
REST_IXS_MODE = "all";  % all/max
% ---

assert(all(ws>0));
t = (0:length(v)-1)/Fs;

[WT, WT_frqs, coi] = cwt(v,Fs,VoicesPerOctave=48,TimeBandwidth=60);
margin = max([find(coi<COI_THRESH,1,'first'), ceil((2*pi)/min(ws)*Fs)]);

% --------------------------------------------------------------------
% Find rest indices:

rest_mask = true;
for frq = ws(:)'/(2*pi)
    [~, ix] = min(abs(WT_frqs - frq));
    dA = abs(diff(abs(WT(ix,:)))./abs(WT(ix,1:end-1)));
    rest_mask = rest_mask & (dA<dA_THRESH);
end

rest_mask(1:margin) = false;
rest_mask(end-margin:end) = false;

cc = bwconncomp(rest_mask,4);
rest_ixs = false(size(t));
switch REST_IXS_MODE
    case "max"
        [~,ccIx] = max(cellfun(@length,cc.PixelIdxList));
        if length(cc.PixelIdxList{ccIx})>MIN_REST_DUR*Fs
            rest_ixs(cc.PixelIdxList{ccIx}) = true;
        end
    case "all"
        for ccIx = 1 : length(cc.PixelIdxList)
            if length(cc.PixelIdxList{ccIx})>MIN_REST_DUR*Fs
                rest_ixs(cc.PixelIdxList{ccIx}) = true;
            end
        end
    otherwise
        error("Unknown rest windows mode");
end
assert(any(rest_ixs));


% --------------------------------------------------------------------
% Calc:

ws_refined = zeros(size(ws));
ph = zeros(size(ws));

WT_amp = abs(WT);
WT_phase = unwrap(atan2(imag(WT), real(WT)),pi,2) + pi/2;

assert(~any(isnan(WT_amp(:))));
assert(~any(isnan(WT_phase(:))));

for i = 1 : length(ws)

    [~, ix] = min(abs(WT_frqs - ws(i)/(2*pi)));

    % --- Refine omega:
    % Method: 
    % Assuming a constant phase offset (=phi), at rest the WT phase = w*t + phi + noise/
    % Therefore, fitting a line to WT phase, will give us w & phi.
    % The line is fitted using ransac. The ransac-margin is taken to be the
    % MAD relative to a naive estimate of w.

    % naive estimate = median of d(phase)/dt:
    w_est = median(diff(WT_phase(ix,rest_ixs)))*Fs;
    % RANSAC, using MAD from naive estimate as margin:
    phase_diff = abs(WT_phase(ix,rest_ixs) - w_est*t(rest_ixs));
    phase_diff_MAD = median(abs(phase_diff-median(phase_diff)));
    ransac_result = ransac_line([t(rest_ixs)',WT_phase(ix,rest_ixs)'], phase_diff_MAD);

    assert(abs(1-ransac_result(1)/ws(i)) < W_REFINE_TOL);  % sanity

    ws_refined(i) = ransac_result(1);   % omega
    ph(i) = ransac_result(2);           % base angular offset

end
