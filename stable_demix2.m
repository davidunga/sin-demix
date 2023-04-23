function dmx = do_demix(v, Fs, ws, options)

arguments
    v
    Fs
    ws
    options.demix_mode = "wi"
    options.stable = true
    options.demix_params = {}
end

demix_fnc = str2func(string(options.demix_mode) + "_demix");

% demix as usual:
dmx0 = demix_fnc(v, Fs, ws, options.demix_params{:});

if ~options.stable
    dmx = dmx0;
    return;
end

% susbtract DC, demix again:
dmx = demix_fnc(v(:)' - dmx0.dc, Fs, ws, options.demix_params{:});

% final result: treat DC component of "no DC" demix as error, and therefore substract 
% if from the DC component. Take sin components from "no DC" demix as-is:
dmx.dc = dmx0.dc - dmx.dc;
