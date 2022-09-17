function r = randrng(rng, varargin)
% random within range

assert(length(rng) == 2);
assert(rng(2) > rng(1));
r = diff(rng) * rand(varargin{:}) + rng(1);
