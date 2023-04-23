function ii = randstep(sz,n,seed)

if length(sz)==1, sz=[1,sz]; end
if nargin==2, seed=1; end
rng(seed);

ii = zeros(sz);
ixs = sort(randperm(length(ii), 2 * n));
for i = 1:2:(2*n)
    ifm = ixs(i);
    ito = ixs(i+1);
    ii(ifm:ito) = 1;
end