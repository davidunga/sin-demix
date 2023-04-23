function dmx = do_demix(v, Fs, ws, opts)

arguments
    v
    Fs
    ws
    opts.mode = "wi"
    opts.stable = true
    opts.params = {}
end

demix_fnc = str2func(string(opts.mode) + "_demix");

% demix as usual:
dmx0 = demix_fnc(v, Fs, ws, opts.params{:});

if ~opts.stable
    dmx = dmx0;
    return;
end

% susbtract DC, demix again:
dmx = demix_fnc(v(:)' - dmx0.dc, Fs, ws, opts.params{:});

% final result: treat DC component of "no DC" demix as error, and therefore substract 
% if from the DC component. Take sin components from "no DC" demix as-is:
dmx.dc = dmx0.dc - dmx.dc;
