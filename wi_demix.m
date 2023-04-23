function dmx = wi_demix(v,Fs,ws)

arguments
    v
    Fs
    ws
end

assert(all(ws>0));
t = (0:length(v)-1)/Fs;
a = nan([length(ws), length(t)]);
p = nan(size(a));

frq_lims = [.5*min(ws), 2*max(ws)]/(2*pi);

TimeBandwidth_for_w1 = 60;

for i = 1:length(ws)
    tbw = TimeBandwidth_for_w1 * (ws(i) / ws(1));
[WT, WT_frqs] = cwt(v,Fs,VoicesPerOctave=48, TimeBandwidth=tbw, FrequencyLimits=frq_lims);

WT_amp = abs(WT);
WT_phase = unwrap(atan2(imag(WT), real(WT)),pi,2) + pi/2;
[WT_tGrid, WT_fGrid] = meshgrid(t, WT_frqs);


%for i = 1 : length(ws)
    [tt, ff] = meshgrid(t, ws(i)/(2*pi));
    a(i,:) = interp2(WT_tGrid, WT_fGrid, WT_amp, tt, ff, "spline");
    p(i,:) = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline") - ws(i)*t;
end

dmx = MixParams(Fs=Fs, ws=ws, v=v, a=a, p=p);
1;
% 
% 
% [pw1,ii1]=peak_width(sum(WT_amp,2),ws(1)/(2*pi),WT_frqs);
% [pw2,ii2]=peak_width(sum(WT_amp,2),ws(2)/(2*pi),WT_frqs);
% 
% nfs = round(2*pw2);
% [tt1, ff1] = meshgrid(t, linspace(WT_frqs(ii1(end)), WT_frqs(ii1(1)), nfs));
% [tt2, ff2] = meshgrid(t, linspace(WT_frqs(ii2(end)), WT_frqs(ii2(1)), nfs));
% a1 = interp2(WT_tGrid, WT_fGrid, WT_amp, tt1, ff1, "spline");
% a2 = interp2(WT_tGrid, WT_fGrid, WT_amp, tt2, ff2, "spline");
% 
% aa1 = movavg(a1, 8*2*pi/ws(2)*Fs);
% aa2 = movavg(a2, 8*2*pi/ws(2)*Fs);
% 
% [~,pth1]=max(aa1,[],1);
% [~,pth2]=max(aa2,[],1);
% 
% inds1=sub2ind(size(a1),pth1,1:length(pth1));
% inds2=sub2ind(size(a2),pth2,1:length(pth2));
% plt a1(inds1) a2(inds2)
% 
% a = nan([length(ws), length(t)]);
% p = nan(size(a));
% 
% fm = mean(ws)/(2*pi);
% nfs = 2*sum(WT_frqs>fm);
% [tt1, ff1] = meshgrid(t, linspace(min(WT_frqs), fm, nfs));
% [tt2, ff2] = meshgrid(t, linspace(fm, max(WT_frqs), nfs));
% 
% a1 = interp2(WT_tGrid, WT_fGrid, WT_amp, tt1, ff1, "spline");
% a2 = interp2(WT_tGrid, WT_fGrid, WT_amp, tt2, ff2, "spline");
% 
% for i = 1 : length(ws)
%     [tt, ff] = meshgrid(t, ws(i)/(2*pi));
%     a(i,:) = interp2(WT_tGrid, WT_fGrid, WT_amp, tt, ff, "spline");
%     p(i,:) = interp2(WT_tGrid, WT_fGrid, WT_phase, tt, ff, "spline") - ws(i)*t;
% end



