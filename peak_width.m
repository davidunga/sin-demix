function [w,ixs] = peak_width(y,px,x,options)

arguments
    y
    px
    x = []
    options.sd = false
end

if isempty(x)
    x = 1 : length(y);
end
[~,ix] = min(abs(x-px));

factor = .5;
ix1 = find(y(1:ix)<factor*y(ix),1,'last');
ix2 = find(y(ix:end)<factor*y(ix),1,'first') + ix - 1;
w = abs(x(ix2)-x(ix1));
if options.sd
    w = w / 2.355;
end
ixs = [ix1,ix,ix2];
