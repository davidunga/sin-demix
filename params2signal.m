function [v, comps] = params2signal(params, t)
% build signal from params struct
% INPUT:
%   params - params struct
%   t - time vector
% OUTPUT:
%   v - the mixed signal, size = [1, length(t)]
%   comps - components (dc, sin1, sin2), size = [3, length(t)]

L = length(t);
a0 = expand(params.a0, L);
a1 = expand(params.a1, L);
a2 = expand(params.a2, L);
w1 = expand(params.w1, L);
w2 = expand(params.w2, L);
p1 = expand(params.p1, L);
p2 = expand(params.p2, L);

comps = [a0; a1 .* sin(w1 .* t + p1); a2 .* sin(w2 .* t + p2)];
v = sum(comps, 1);

function x = expand(x, L)
if length(x) == 1
    x = x * ones(1,L);
end
assert(length(x) == L);
