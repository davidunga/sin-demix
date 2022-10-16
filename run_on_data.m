function run_on_data(DataBase, options)

% ---
% handle input:

arguments
    DataBase            % either struct or matfile
    options.ds = true   % flag - downsample to nyquist?
end

if ~isstruct(DataBase)
    assert(isstring(DataBase) || ischar(DataBase), "DataBase should be a struct or matfile");
    db = load(DataBase);
    DataBase = db.DataBase;
    clear db;
end

% ---
% parse data:

w1 = DataBase.FR.ff * 2 * pi;
w2 = DataBase.FR.ff2 * 2 * pi;
v = DataBase.FR.V(:);
Fs = DataBase.FR.sr;
t = (0:(length(v)-1))'/Fs;

if options.ds
    % downsample
    Fs = 2 * max(w1,w2);
    v = interp1(t, v, linspace(t(1), t(end), round((t(end) - t(1)) * Fs)));
    t = linspace(t(1), t(end), round((t(end) - t(1)) * Fs));
end

% ---
% demix:

comps_hat = demix(v, t, w1, w2);

% ---
% plot:

draw_components(t,comps_hat,v);
sgtitle(DataBase.FR.FileName,"Interpreter","none");
