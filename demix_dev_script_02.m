% demix_dev_script_02

clear;
METHOD = "corr";

Fs = 1000;
ws = [100, 200];

dur = 5;
N = round(dur * Fs);

dc = piecewise_const(N, [0, 2, .5], [0, .25, .75]);
amp2 = piecewise_const(N, [5, 4.5], [0, .4]);
amp1 = piecewise_const(N, [6.5, 7.0], [0, .6]);

gt = MixParams(Fs=Fs,ws=ws,dc=dc,a={amp1,amp2},ph=[45, 25]*pi/180);
t = gt.t;

% -----

[ws,~,margin] = refine_sine_params(gt.v, Fs, ws);
if METHOD == "corr"
    dmx = wt_demix2c(gt.v, Fs, ws);
elseif METHOD == "linear"
    dmx = linear_demix2(gt.v, Fs, ws);
else
    error("Unknown method");
end

GT_COLORS = ["r","r","m"];
EST_COLORS = ["k","b","c"];

figure();
tl = tiledlayout(2,1,TileSpacing="compact",Padding="compact");
title(tl,sprintf("w1=%3.2f w2=%3.2f Fs=%3.2f method=%s", ws(1), ws(2), Fs, METHOD));

gt_amps = gt.amps;
dmx_amps = dmx.amps;
nexttile(); hold on;
for i = 1 : size(gt_amps,1)
    plot(t, gt_amps(i,:), GT_COLORS(i), DisplayName=sprintf("a%d TRUE",i-1));
    plot(t, dmx_amps(i,:), EST_COLORS(i), DisplayName=sprintf("a%d EST",i-1));
end
grid minor;
legend();
xlabel("t [sec]");
title("Amplitudes");

nexttile(); hold on;
for i = 1 : size(gt.p,1)
    plot(t, gt.p(i,:), GT_COLORS(i+1), DisplayName=sprintf("p%d TRUE",i));
    plot(t, dmx.p(i,:), EST_COLORS(i+1), DisplayName=sprintf("p%d EST",i));
end
grid minor;
legend();
xlabel("t [sec]");
title("Phases");

linkaxes(allaxes(),'x');
xlim(t([margin,end-margin]));



figure();
tiledlayout(2,1,TileSpacing="compact",Padding="compact");

c = dmx.comps;
nexttile(); hold on;
plot(gt.t,gt.v-sum(c(2:end,:),1),'r',DisplayName='DC by residual');
plot(gt.t,dmx.dc,'b',DisplayName='Estimated DC');
grid minor;
legend();
xlabel("t [sec]");
title("DCs");

nexttile(); hold on;
plot(gt.t,dmx.dc-(gt.v-sum(c(2:end,:),1)),'k');
grid minor;
legend();
xlabel("t [sec]");
title("EstDC - ResidualDC");

linkaxes(allaxes(),'x');
xlim(t([margin,end-margin]));