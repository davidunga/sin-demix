function [comps,opt] = wt_demix4(v,Fs,ws)

arguments
    v
    Fs
    ws
end

assert(length(ws)==1 || all(diff(ws) > 0), "Frequencies must be monotonically increasing");

assert(ws(1)==0);
w1 = ws(2);
w2 = ws(3);

params = struct();
%params.n = 5;
params.dur = 2*pi/min([w1,w2]);
params = namedargs2cell(params);

[wv1,~] = make_naive_wavelet(Fs,w1/(2*pi),params{:});
[wv2,~] = make_naive_wavelet(Fs,w2/(2*pi),params{:});
t = 0:(1/Fs):2*2*pi/min([w1,w2]);

w11 = max(conv(wv1,sin(w1*t),"same"));
w12 = max(conv(wv1/w11,sin(w2*t),"same"));
w22 = max(conv(wv2,sin(w2*t),"same"));
w21 = max(conv(wv2/w22,sin(w1*t),"same"));

s1 = conv(v,wv1/w11,"same");
s2 = conv(v,wv2/w22,"same");

a1 = (amplitude(s1) - w12);
a2 = (amplitude(s2) - w21);

c1 = set_amplitude(s1,a1);
c2 = set_amplitude(s2,a2);

comps = [v-(s1+s2);s1;s2];

opt.r_dur = 2*pi/max(ws);
opt.r = round(opt.r_dur * Fs);
