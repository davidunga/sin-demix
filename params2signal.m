function [v, comps] = params2signal(params, t, smooth, r)
% build signal from params struct
% INPUT:
%   params - params struct
%   t - time vector
%   smooth - logical, apply smoothing? default = false
%   r - margin - replicate first and last [r] elements, default = 0
% OUTPUT:
%   v - the mixed signal, size = [1, length(t)]
%   comps - components (dc, sin1, sin2), size = [3, length(t)]

if ~exist('smooth', 'var'), smooth=false; end
if ~exist('r', 'var'), r=0; end

t = t(:)';
L = length(t);
a0 = expand(params.a0, L, r);
a1 = expand(params.a1, L, r);
a2 = expand(params.a2, L, r);
w1 = expand(params.w1, L, r);
w2 = expand(params.w2, L, r);
p1 = expand(params.p1, L, r);
p2 = expand(params.p2, L, r);

if smooth
    Fs = (L-1)/(t(end) - t(1));
    win = 1 * (2*pi)/max(w1(1),w2(1)) * Fs;
    a0 = smoothdata(a0,"movmean",win);
end

comps = [a0; a1 .* sin(w1 .* t + p1); a2 .* sin(w2 .* t + p2)];

if r > 0
    comps(:,1:r) = repmat(comps(:,r+1),[1,r]);
    comps(:,end-(r-1):end) = repmat(comps(:,end-r),[1,r]);
end

v = sum(comps, 1);

function x = expand(x, L, r)
if length(x) == 1
    x = x * ones(1,L);
elseif r > 0
    x(1:r) = x(r+1);
    x(end-r:end) = x(end-r);
end
assert(length(x) == L);
