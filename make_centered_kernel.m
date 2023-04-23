function kernel = make_centered_kernel(k)
% Given a kernel size [k], where k can be fractional, returns a symmetric odd-length
% kernel of the form: [p,..,1,1,1,..,p], where 0<p<1, s.t. sum(kernel)==k
% This kernel can be used to simulate a kernel of size [k], while
% maintaining the properties of an odd-length kernel.
% examples:
%   make_centered_kernel(2) -> [.5,1,.5]
%   make_centered_kernel(2.5) -> [.75,1,.75]
%   make_centered_kernel(3) -> [1,1,1]

assert(k>0);

if k < 1
    kernel = k;
    return;
end

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

assert(abs(k-sum(kernel)) < 10*eps*k);
