function ret = calc_signal_properties_CWT(v, Fs, opts)
% Given a sine-mixture signal [v] sampled at [Fs], returns a struct with fields:
%   - ws - Estimate of sine omegas
%   - ph - Estimate of base phase offset for each omega
%   - margin - margin of influence [index count]
%   - rest_ixs - mask vector, same size as [v], where no amplitude changes occure
%   - cwtRes - result of cwt

arguments
    v
    Fs
    opts.rest_mode = "all"   % "all"/"max" - take all rest windows, or only the maximal duration window
    opts.cwt_params = {'VoicesPerOctave', 48, 'TimeBandwidth', 60, 'FrequencyLimits', [0, 2*pi*1000]}
end

% --------------------------------------------------------------------
% Params

dA_THRESH = .00005; % rest is defined when relative amplitude change is lower than this
MIN_REST_DUR = .05; % exclude rest windows shorter than this [sec]
AMP_RATIO_EPS = .1; % When computing omegas, discard amplitude peaks lower than [AMP_RATIO_EPS]*max_amplitude
AMP_RATIO_LOW = .7; % When computing omegas, assert that the ratio between the min and max peak is at least this value 
FRQ_SEP = 10;       % Minimal separtion between frequencies [Hz]

% --------------------------------------------------------------------
% Find cwt-based estimate of omegas

% compute wavelet transform:
cwtRes = struct();
[cwtRes.WT, cwtRes.frqs, cwtRes.coi] = cwt(v,Fs,opts.cwt_params{:});

% get indices of frequnecy peaks:
amp = sum(abs(cwtRes.WT),2);
frq_ixs = length(amp) + 1 - find(islocalmax(amp(end:-1:1),MinSeparation=FRQ_SEP,SamplePoints=cwtRes.frqs(end:-1:1)));
frq_ixs = frq_ixs(amp(frq_ixs) > AMP_RATIO_EPS*max(amp(frq_ixs)));

assert(min(amp(frq_ixs)) > AMP_RATIO_LOW*max(amp(frq_ixs))); % make sure all peaks are roughly the same magnitude

% influence margin is the COI  that includes all omegas (+tolerance=FRQ_SEP):
margin = find(cwtRes.coi<min(cwtRes.frqs(frq_ixs))-FRQ_SEP,1,'first');


% --------------------------------------------------------------------
% Find rest indices:

% all rest indices
rest_mask = true;
for frq_ix = frq_ixs(:)'
    dA = abs(diff(abs(cwtRes.WT(frq_ix,:)))./abs(cwtRes.WT(frq_ix,1:end-1)));
    rest_mask = rest_mask & (dA<dA_THRESH);
end

rest_mask(1:margin) = false;
rest_mask(end-margin+1:end) = false;

cc = bwconncomp(rest_mask,4); % group rest durations into rest-windows
rest_ixs = false(size(v));
switch opts.rest_mode
    case "max"
        % take only max duration window (if its sufficiently long)
        [~,ccIx] = max(cellfun(@length,cc.PixelIdxList));
        if length(cc.PixelIdxList{ccIx})>MIN_REST_DUR*Fs
            rest_ixs(cc.PixelIdxList{ccIx}) = true;
        end
    case "all"
        % take all sufficiently long windows
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
% Calc refined omegas and phase offset

% Method: 
% Assuming constant phase offset (phi) at rest, the phase is: w*t + phi + noise.
% Therefore, fitting a line to the phase, will give us w & phi.
% The line is fitted using ransac. The ransac-margin is taken to be the
% MAD relative to a initial omegas estimate.

ws = zeros([1,length(frq_ixs)]);
ph = zeros([1,length(frq_ixs)]);

WT_phase = unwrap(atan2(imag(cwtRes.WT), real(cwtRes.WT)),pi,2) + pi/2;
t = (0:length(v)-1)/Fs;

for i = 1:length(frq_ixs)

    frq_ix = frq_ixs(i);
    
    % naive estimate = median of d(phase)/dt:
    w_est = median(diff(WT_phase(frq_ix,rest_ixs)))*Fs;

    % RANSAC, using MAD from naive estimate as margin:
    phase_diff = abs(WT_phase(frq_ix,rest_ixs) - w_est*t(rest_ixs));
    phase_diff_MAD = median(abs(phase_diff-median(phase_diff)));
    ransac_result = ransac_line([t(rest_ixs)',WT_phase(frq_ix,rest_ixs)'], phase_diff_MAD);

    ws(i) = ransac_result(1);   % omega
    ph(i) = ransac_result(2);   % base angular offset

end

ret = vars2struct(["ws", "ph", "margin", "rest_ixs", "cwtRes"]);


% --------------------------------------------------------------------

DBG_PLOT = false;
if DBG_PLOT

    figure();
    tiledlayout(2,1,TileSpacing="compact",Padding="compact");
    nexttile(); hold on;
    plot(cwtRes.frqs, amp, DisplayName="Power");
    plot(cwtRes.frqs(frq_ixs), amp(frq_ixs), 'rs', DisplayName="Omegas from CWT");
    plot(ws/(2*pi), amp(frq_ixs), 'r*', DisplayName="Refined Omegas");
    xlabel("Frq [Hz]");
    ylabel("Power");
    legend();

    nexttile(); hold on;
    avg_amp = mean(abs(cwtRes.WT(frq_ixs,:)),1);
    plot(t, avg_amp, 'b-', DisplayName="Avrg amplitude");
    plot(t(rest_ixs), avg_amp(rest_ixs), 'r.', DisplayName="Rest");
    xlabel("Time [sec]");
    ylabel("Avrg Amplitude");
    legend();

    uiwait();
end

