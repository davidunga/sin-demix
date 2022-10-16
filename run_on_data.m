function run_on_data(DataBase, options)

% --------------------------------------
% handle input:

arguments
    DataBase            % either data struct or path to datafile
    options.ds = true   % flag - downsample to nyquist?
end

if ~isstruct(DataBase)
    load(DataBase,'DataBase');
end

% --------------------------------------
% parse data:

w1 = DataBase.FR.ff * 2 * pi;
w2 = DataBase.FR.ff2 * 2 * pi;
v = DataBase.FR.V(:);

if options.ds
    % downsample
    Fs = round(2 * max(w1,w2));
    v = rsmpl(v,DataBase.FR.sr,Fs);
else
    Fs = DataBase.FR.sr;
end

% --------------------------------------
% demix:

comps_hat = demix(v, Fs, w1, w2);

% --------------------------------------
% plot:

draw_components(Fs,comps_hat,v);
sgtitle(DataBase.FR.FileName,"Interpreter","none");
