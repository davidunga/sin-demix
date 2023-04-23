function [c,dc] = local_sine_corr(x, Fs, w, ph)

d = round(Fs*2*pi/w);

pd = ceil(length(x)/d)*d - length(x);
x = padarray(x(:),pd,nan,'post')';

t = (0:(length(x)-1))/Fs;
s = sin(w*t + ph);
s(end-pd+1:end) = nan;

x = reshape(x, d, [])';
dc = nanmean(x,2);
s = reshape(s, d, [])';

c = nansum((x-dc).*s,2);
c = c ./ nansum(s.^2,2);

c = movmedian(c,3);
dc = movmedian(dc,3);

dc = reshape(repmat(dc,[1,d])',1,[]);
c = reshape(repmat(c,[1,d])',1,[]);

c = c(1:end-pd);
dc = dc(1:end-pd);