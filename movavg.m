function [ret, kernel] = movavg(x,k,options)
% robust moving average using centered kernel.
% supports non-even (odd and fractional) window sizes (k).

arguments
    x
    k
    options.sum = false  % if true, returns moving sum instead of average
    options.dim = 1;
end

kernel = make_centered_kernel(k);
if isvector(x)
    ret = conv(x,kernel,'same');
else
    ret = conv2(x,kernel,'same');
end
if ~options.sum
    ret = ret / sum(kernel);
end

n = length(kernel);
ret(1:floor(n/2)) = ret(ceil(n/2));
ret(end-floor(n/2)+1:end) = ret(end-ceil(n/2)+1);
