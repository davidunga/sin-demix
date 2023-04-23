function [comps, ret] = sine_corr_fixed_phase(v, Fs, ws, ph)

v = v(:)';
d = round(Fs*2*pi/min(ws));
dc = movmean(v,d);
vv = v - dc;

t = (0:(length(v)-1))/Fs;

a = nan([length(ws), length(v)]);
for i = 1 : length(ws)
    s = sin(ws(i).*t + ph(i));
    a(i,:) = movmean(vv.*s,d);
    a(i,:) = a(i,:) ./ movmean(s.*s,d);
end

ret = MixParams(Fs=Fs, ws=ws, dc=dc, a=a, ph=ph);
comps = ret.comps;

return;