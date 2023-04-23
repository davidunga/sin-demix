function [a,u,l,s] = amplitude(x,n)
% amplitude and envelope

if nargin==1, n=1; end

[u,l]=envelope(x,1);


a = .5 * (u - l);
s = .5 * (u + l);
