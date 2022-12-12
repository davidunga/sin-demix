function [comps,opt] = naive_demix(v, Fs, ws)

f1 = ws(2)/(2*pi);
f2 = ws(3)/(2*pi);

T = min(f1,f2)^-1;
[wv1,~] = make_naive_wavelet(Fs,f1,T);
[wv2,~] = make_naive_wavelet(Fs,f2,T);

mag1 = max(conv(sin(2*pi*f1*(0:4*T*Fs)/Fs),wv1,"same"));
mag2 = max(conv(sin(2*pi*f2*(0:4*T*Fs)/Fs),wv2,"same"));

wv1 = wv1 / mag1;
wv2 = wv2 / mag2;

c1 = conv(v,wv1,"same");
c2 = conv(v,wv2,"same");
comps = [v-(c1+c2);c1;c2];

opt.r_dur = T;
opt.r = round(T*Fs);