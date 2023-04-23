function dmx = wi_demix_DEP(v,Fs,ws,opts)

arguments
    v
    Fs
    ws
    opts.cwtRes = []
    opts.cwt_params = {'VoicesPerOctave', 48, 'TimeBandwidth', 60, 'FrequencyLimits', [0, 2*pi*1000]}
end

assert(all(ws>0));
t = (0:length(v)-1)/Fs;

% --- Calc:

a = nan([length(ws), length(t)]);
p = nan(size(a));

if ~isempty(opts.cwtRes)
    cwtRes = opts.cwtRes;
    clear options.cwtRes;
else
    cwtRes = struct();
    [cwtRes.WT, cwtRes.frqs] = cwt(v,Fs,opts.cwt_params{:});
end

WT_amp = abs(cwtRes.WT);
WT_phase = unwrap(atan2(imag(cwtRes.WT), real(cwtRes.WT)),pi,2) + pi/2;
[WT_tGrid, WT_fGrid] = meshgrid(t, cwtRes.frqs);

assert(~any(isnan(WT_amp(:))));
assert(~any(isnan(WT_phase(:))));

for i = 1 : length(ws)
    % Get amplitude and phase at refined omega:
    [tt, ff] = meshgrid(t, ws(i)/(2*pi));
    a(i,:) = interp2(WT_tGrid, WT_fGrid, WT_amp, tt, ff, "spline");
    p(i,:) = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline") - ws(i)*t;
end

dmx = MixParams(Fs=Fs, ws=ws, v=v, a=a, p=p);
