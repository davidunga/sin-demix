function [comps,opt] = wt_demix(v,Fs,ws,options)

arguments
    v
    Fs
    ws
    options.params = [2,50];
    options.normP2 = true;
end

assert(length(ws)==1 || all(diff(ws) > 0), "Frequencies must be monotonically increasing");

has_dc = ws(1) == 0;
ws = ws(ws>0);

WF = frqlocwt(v,Fs,ws/(2*pi),options.params,options.normP2);
comps = real(WF);

if has_dc
    comps = [v - sum(comps,1); comps];
end

opt.r_dur = 2*pi/max(ws);
opt.r = round(opt.r_dur * Fs);
