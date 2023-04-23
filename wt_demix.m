function dmx = wt_demix(v,Fs,ws,opts)

arguments
    v
    Fs
    ws
    opts.params = [3,40];
end

assert(all(ws>0));

if opts.params(2) == -1
    opts.params(2) = max(morseBandWidth_Time2Power(ws/(2*pi), opts.params(1) * 39.99));
    opts.params
end
opts.params = "morse";
wts = cwt_tight(v,Fs,ws/(2*pi),wv=opts.params,r=3);

t = (0:length(v)-1)/Fs;
a = nan([length(ws), length(t)]);
p = nan(size(a));
for i = 1 : length(ws)
    WT = wts(i).WT;
    WT_frqs = wts(i).frqs;
    WT_amp = abs(WT);
    WT_phase = unwrap(atan2(imag(WT), real(WT)),pi,2) + pi/2;
    [WT_tGrid, WT_fGrid] = meshgrid(t, WT_frqs);
    [tt, ff] = meshgrid(t, ws(i)/(2*pi));
%     [tt, ff] = meshgrid(t, ws(i)/(2*pi)+linspace(-3,3,500));
%     aa = interp2(WT_tGrid, WT_fGrid, WT_amp, tt, ff, "spline");
%     pp = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline");
%     [mx,mxi] = max(aa,[],1);
%     mxi = round(movavg(mxi,2*pi/ws(1)*Fs));
%     a(i,:) = aa(sub2ind(size(pp),mxi,1:length(mxi)));
%     p(i,:) = pp(sub2ind(size(pp),mxi,1:length(mxi))) - ws(i)*t;
    a(i,:) = interp2(WT_tGrid, WT_fGrid, WT_amp, tt, ff, "spline");
    p(i,:) = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline") - ws(i)*t;
end

dmx = MixParams(Fs=Fs, ws=ws, v=v, a=a, p=p);

