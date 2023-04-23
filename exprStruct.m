function s = exprStruct(params)
% Make struct of experiment parameters

arguments
    params.Fs               % Sampling rate [Hz]
    params.ws = []          % Injected omegas [Rad/sec]
    params.rest_win = []    % Rest window start & end times [sec]
end

s = struct(Fs=params.Fs, ws=sort(params.ws), rest_win=params.rest_win);

if ~isempty(s.ws)
    assert(all(s.ws>0));
end

if ~isempty(s.rest_win)
    assert(numel(s.rest_win)==2);
    assert(diff(s.rest_win)>0);
end
