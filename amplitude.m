function [a,u,l] = amplitude(x)
% amplitude and envelope


ii = islocalmax(x);
uu = interp1(find(ii),x(ii),1:length(x),"pchip");
ii = islocalmax(-x);
ll = interp1(find(ii),x(ii),1:length(x),"pchip");
u = max(uu, ll);
l = min(uu, ll);
a = .5 * (u - l);