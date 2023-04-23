function [comps,opt] = adjust_to_fixed_phase(v,Fs,opt)

a = opt.a;
p = opt.p;
ws = opt.ws_refined;
ph = opt.ph;
t = (0:(length(v)-1))/Fs;

[c1,dc1] = local_sine_corr(v, Fs, ws(1), ph(1));
[c2,dc2] = local_sine_corr(v, Fs, ws(2), ph(2));

rest_ixs = (t>.35) & (t<.65);
j = 2;
wss = linspace(.99,1.01,501)*ws(j);
phs = linspace(.98,1.02,101)*ph(j);
sd = nan(size(wss));
best_sd = 100000;
for i = 1 : length(wss)
    for jj = 1 : length(phs)
        [cc,dcc] = local_sine_corr(v, Fs, wss(i), phs(jj));
        sd = std(cc(rest_ixs));
        if sd < best_sd
            best_sd = sd;
            best_i = i;
            best_jj = jj;
        end
    end
end
[~,mni] = min(sd);

comps = zeros([size(a,1)+1,length(v)]);

for i = 1 : length(ph)
    time_dep_phase = p(i,:) - ph(i);
    a(i,:) = a(i,:) .* cos(time_dep_phase);
    p(i,:) = ph(i);
    comps(i+1,:) = a(i,:).*sin(ws(i)*t + p(i,:));
end

comps(1,:) = v - sum(comps(2:end,:),1);
opt.a = a;
opt.p = p;
opt.ph = ph;
