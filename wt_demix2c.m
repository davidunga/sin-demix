function mp = wt_demix2c(v,Fs,ws)

arguments
    v
    Fs
    ws
end

t = (0:length(v)-1)/Fs;
a = nan([length(ws), length(t)]);
p = nan(size(a));

win_sz = Fs*2*pi/min(ws);
dc = movavg(v, win_sz);

for i = 1 : length(ws)
    z = movavg(v.*exp(1i*ws(i)*t),win_sz);
    a(i,:) = 2*abs(z);
    p(i,:) = atan2(real(z), imag(z));
end

mp = MixParams(Fs=Fs, ws=ws, dc=dc, a=a, p=p);
