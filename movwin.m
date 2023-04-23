function ret = movwin(x,k,options)
% robust moving window.
% window size (k) is adjusted to be an odd integer, and the filtering
% weights are computed to simulate [k].

arguments
    x
    k
    options.mode = "mean"
end

assert(any(strcmp(options.mode, ["sum", "mean"])));

if round(k)==k
    if mod(k,2) == 0
        n = k + 1;
    else
        n = k;
    end
else
    n = ceil(k);
    if mod(n, 2) == 0
        n = n + 1;
    end
end

kernel = ones([1,n]);
kernel(1) = 1-(n-k)/2;
kernel(end) = 1-(n-k)/2;

assert(sum(kernel)==k);

ret = conv(x,kernel,'same');
if strcmp(options.mode, "mean")
    ret = ret / sum(kernel);
end

ret(1:floor(n/2)) = ret(ceil(n/2));
ret(end-floor(n/2)+1:end) = ret(end-ceil(n/2)+1);