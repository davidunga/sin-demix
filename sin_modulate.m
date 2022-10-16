function result = sin_modulate(t, baseval, w, amp)
% Create an oscillatory time sequence by sin-modulating a base value
% EXAMPLE:
%   Create 10% oscillations around 3, with angular frequency 60 Hz:
%       sin_modulate(3, t, 60, .1);

result = baseval * (1 + amp * sin(w*t));