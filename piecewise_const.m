function x = piecewise_const(L, vals, locs)
% Make a piecewise constant vector.
% INPUT:
%   L - length of vector
%   vals - array of values for each piece.
%   locs - location of piece starts. can be, either:
%       1. An array of indices (between 1 and L), same size as [vals]. First value must be 1.
%       2. An array of relative indices (between 0 and 1), same size as [vals]. First value must be 0.
%       3. 'u' - uniformly spaced pieces (default).
%       4. 'r' - randomly spaced pieces.
% OUTPUT:
%   array [1,L]

if ~exist('locs', 'var'), locs = 'u'; end

if ~isnumeric(locs)
    if locs == 'u'
        locs = round(linspace(1, L, length(vals) + 1));
    elseif locs == 'r'
        locs = [1, 1 + sort(randperm(L - 2, length(vals) - 1)), L];
    else
        error('Unknown locs value');
    end
else
    if locs(1) == 0
        assert(all(locs <= 1));
        locs = round(locs * (L - 1)) + 1;
    end
    assert(locs(1) == 1, "First location must be 1");
    locs(end+1) = L;
end


x = zeros([1, L]);
for i = 1 : length(vals)
    x(locs(i):locs(i+1)) = vals(i);
end
