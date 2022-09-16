function result = sin_modulate(baseval, t, w, amp)
% Create an oscillatory time sequence by sin-modulating a base value
% INPUT:
%   baseval - scalar
%   t - time vec
%   w - angular frequency. if given as a range, a value is randomly
%       selected from that range
%   amp - amplitude relative to base value. if given as a range, a value is randomly
%       selected from that range.
% OUTPUT:
%   vector, same size as [t]

if length(w) == 2
    w = rand() * diff(w) + w(1);
end
if length(amp) == 2
    amp = rand() * diff(amp) + amp(1);
end

modulator = 1 + amp * sin(w*t/(2*pi));
result = baseval * modulator;