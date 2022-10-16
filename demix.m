function comps_hat = demix(v, t, w1, w2)
% Demix time sequence into a time-modulated sum of sins.
% Specifically, for a signal v(t), assumed to be a mixture of the form:
%   v(t) = a0(t) + a1(t)*sin(w1*t + p1(t)) + a2(t)*sin(w2*t + p2(t))
% Finds [a0(t)], [a1(t)*sin(w1*t + p1(t))], and [a2(t)*sin(w2*t + p2(t))]
% under the assumption that they are varying slower than max(w1,w2).
% INPUT:
%   v, t - signal and time vectors to demix
%   w1, w2 - frequency of the underlying sins
% OUTPUT:
%   estimated components, given as a [3, length(t)] array.

v = v(:);
t = t(:);

A = [ones(size(t)), sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t)];

x = nan(size(A));

Fs = (length(t)-1)/(t(end)-t(1));
win_dur = 2*pi/max(w1,w2);
r = round(win_dur * Fs);
r = r + (1-mod(r,2));
max_eqs = 5;
sample_ixs = round(linspace(-r,r,min(2*r+1,max_eqs)));

for i = (r+1):(size(A,1)-r)
    x(i, :) = A(i + sample_ixs, :) \ v(i + sample_ixs);
end

x(1:r,:) = repmat(x(r+1,:),[r,1]);
x(end-(r-1):end,:) = repmat(x(end-r,:),[r,1]);
x = x';

params_hat = struct();
params_hat.a0 = x(1, :);
params_hat.a1 = sqrt(x(2,:).^2 + x(3,:).^2);
params_hat.p1 = atan2(x(3,:), x(2,:));
params_hat.a2 = sqrt(x(4,:).^2 + x(5,:).^2);
params_hat.p2 = atan2(x(5,:), x(4,:));
params_hat.w1 = w1;
params_hat.w2 = w2;

[~, comps_hat] = params2signal(params_hat, t, false);
