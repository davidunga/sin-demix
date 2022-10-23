function x = rsmpl(x, Fs, newFs, options)
% Resample x from sampling rate [Fs] to [newFs]

arguments
    x
    Fs
    newFs
    options.method = "linear"
end

N = length(x)*newFs/Fs;
if newFs > Fs
    N = ceil(N);
else
    N = floor(N);
end
t_new = (0:(N-1)) / newFs;

t = (0:(length(x)-1)) / Fs;
x = interp1(t,x,t_new,options.method);