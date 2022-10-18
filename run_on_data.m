function result = run_on_data(DataBase, src, options)

% --------------------------------------
% handle input:

arguments
    DataBase                % either data struct or path to datafile
    src = "V"               % "V"/"I"
    options.plot = true     % plot result?
    options.ds = false      % downsample to nyquist?
end

if ~isstruct(DataBase)
    load(DataBase,'DataBase');
end

% --------------------------------------
% parse data:

w1 = DataBase.FR.ff * 2 * pi;
w2 = DataBase.FR.ff2 * 2 * pi;
v = reshape(DataBase.FR.(src),[],1);

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

result = struct();
result.src = src;
result.Fs = Fs;
result.dc = comps_hat(1,:);
result.sin1 = comps_hat(2,:);
result.sin2 = comps_hat(3,:);

if options.plot
    draw_components(result,v);
    sgtitle(DataBase.FR.FileName,"Interpreter","none");
end
