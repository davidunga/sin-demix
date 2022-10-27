function [a,u,l] = amplitude(x,n)
% amplitude and envelope


ii = find(islocalmax(x));
ii = ii(1:n:end);
uu = interp1(ii,x(ii),1:length(x),"pchip");

ii = find(islocalmax(-x));
ii = ii(1:n:end);

ll = interp1(ii,x(ii),1:length(x),"pchip");
u = max(uu, ll);
l = min(uu, ll);
a = .5 * (u - l);