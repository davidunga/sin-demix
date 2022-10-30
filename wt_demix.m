function [comps,opt] = wt_demix(v,Fs,ws,options)

arguments
    v
    Fs
    ws
    options.params = [3,-1];
end

assert(length(ws)==1 || all(diff(ws) > 0), "Frequencies must be monotonically increasing");

has_dc = ws(1) == 0;
ws = ws(ws>0);

if options.params(2) == -1
    options.params(2) = max(morseBandWidth_Time2Power(ws/(2*pi), options.params(1) * 39.99));
    options.params
end

WF = wtbandpass(v,Fs,ws/(2*pi),options.params);
comps = real(WF);

if has_dc
    comps = [v - sum(comps,1); comps];
end

opt.r_dur = 2*pi/max(ws);
opt.r = round(opt.r_dur * Fs);
