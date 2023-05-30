function dmx = wt_demix(v,Fs,ws,opts)

arguments
    v
    Fs
    ws
    opts.params = [3,-1];
end

assert(all(ws>0));

if opts.params(2) == -1
    opts.params(2) = max(morseBandWidth_Time2Power(ws/(2*pi), opts.params(1) * 39.99));
    opts.params
end
%opts.params = "morse";
wts = cwt_tight(v,Fs,ws/(2*pi),wv=opts.params,r=3);

t = (0:length(v)-1)/Fs;
a = nan([length(ws), length(t)]);
p = nan(size(a));
for i = 1 : length(ws)
    WT_phase = unwrap(atan2(imag(wts(i).WT), real(wts(i).WT)),pi,2) + pi/2;
    [WT_tGrid, WT_fGrid] = meshgrid(t, wts(i).frqs);
    [tt, ff] = meshgrid(t, ws(i)/(2*pi));
    a(i,:) = interp2(WT_tGrid, WT_fGrid, abs(wts(i).WT), tt, ff, "spline");
    p(i,:) = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline") - ws(i)*t;
end

dmx = MixParams(Fs=Fs, ws=ws, v=v, a=a, p=p);

