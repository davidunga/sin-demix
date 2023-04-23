function c = movmean_f(x,k)

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

c = conv(x,k);