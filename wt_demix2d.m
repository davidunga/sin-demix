function dmx = wt_demix2d(v,Fs,ws,options)

arguments
    v
    Fs
    ws
    options.cwtRes = []
end

assert(all(ws>0));
t = (0:length(v)-1)/Fs;

% --- Calc:

a = nan([length(ws), length(t)]);
p = nan(size(a));

if ~isempty(options.cwtRes)
    cwtRes = options.cwtRes;
    clear options.cwtRes;
else
    cwtRes = struct();
    [cwtRes.WT, cwtRes.frqs] = cwt(v,cwtRes.Fs,VoicesPerOctave=48,TimeBandwidth=60);
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
