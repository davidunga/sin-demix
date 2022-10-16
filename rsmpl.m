function x = rsmpl(x, Fs, newFs, options)
% resample x from Fs newFs

arguments
    x
    Fs
    newFs
    options.method = "linear"
end

t = (0:(length(x)-1)) / Fs;
tnew = (0:((round(length(x)*newFs/Fs))-1)) / newFs;
x = interp1(t,x,tnew,options.method);