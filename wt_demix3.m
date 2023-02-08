function [comps,opt] = wt_demix3(v,Fs,ws)

arguments
    v
    Fs
    ws
end

assert(length(ws)==1 || all(diff(ws) > 0), "Frequencies must be monotonically increasing");

has_dc = ws(1) == 0;
ws = ws(ws>0);

WF = [  morlet_tform(v,Fs,ws(1)/(2*pi),n=4);
        morlet_tform(v,Fs,ws(2)/(2*pi),n=4)];

comps = real(WF);

if has_dc
    comps = [v - sum(comps,1); comps];
end

opt.r_dur = 2*pi/max(ws);
opt.r = round(opt.r_dur * Fs);
