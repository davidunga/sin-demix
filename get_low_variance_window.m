function ixs = get_low_variance_window(x, Fs, ws, opts)
% get indices of window of duration=[rest_dur], where the
% std of V0 was minimal.

arguments
    x
    Fs
    ws
    opts.cycles_in_win = 25
    opts.margin_cycles = 5
end

cycle_dur = 2*pi/min(ws);
win_dur = opts.cycles_in_win*cycle_dur;
r = round(win_dur*Fs/2);
margin = 3*r;
xx = movstd(x,2*r+1);
[~,mni] = min(movavg(xx(margin:end-margin), 2*r+1));
ixs = (mni+margin-r):(mni+margin+r);
return