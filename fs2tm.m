function t = fs2tm(Fs,dur,sz)
% time vec from sampling rate and duration
% input: either- fs2tm(Fs,dur) or fs2tm(Fs,[],sz)

if nargin == 3
    assert(isempty(dur));
    if length(sz) == 1
        sz = [sz,1];
    end
else
    sz = [round(Fs*dur),1];
end

t = reshape(0:(prod(sz)-1), sz)/Fs;