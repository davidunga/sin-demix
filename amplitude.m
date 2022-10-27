function [a,u,l,s,ui,li] = amplitude(x,n)
% amplitude and envelope

if nargin==1, n=1; end

[u,ui] = get_upper_bound(x,n);
[l,li] = get_upper_bound(-x,n);
l = -l;

a = .5 * (u - l);
s = .5 * (u + l);

function [y,ii] = get_upper_bound(x,n)
ii = false(size(x));
locmax = find(islocalmax(x));
ii(locmax(1:n:end)) = true;
y = interp1(locmax(1:n:end),x(ii),1:length(x),"pchip");
