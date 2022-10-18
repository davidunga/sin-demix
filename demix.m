function [comps_hat,err] = demix(v, Fs, w1, w2)
% Demix time sequence into a time-modulated sine-like components.
% Specifically, for a signal v(t), assumed to be a mixture of the form:
%   v(t) = a0(t) + a1(t)*sin(w1*t + p1(t)) + a2(t)*sin(w2*t + p2(t)), with
%   known frequencies w1, w2, finds the components: 
%   c1(t) = a0(t)
%   c2(t) = a1(t)*sin(w1*t + p1(t))
%   c3(t) = a2(t)*sin(w2*t + p2(t))
% Under the assumption that they are varying slower than max(w1,w2).
% INPUT:
%   v, Fs - signal to demix and its sampling rate
%   w1, w2 - frequency of the underlying sines
% OUTPUT:
%   comps_hat - estimated components, given as a [3, length(v)] array.
%   err - reconstruction error- MSE normalized by signal variance, excluding boundaries

v = v(:);
t = (0:(length(v)-1))'/Fs;

A = [ones(size(t)), sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t)];

r_dur = 2*pi/max(w1,w2); % window radius, in seconds
r = round(r_dur * Fs); % sampling radius
sample_ixs = round(linspace(-r,r,5)); % 0-centered sampling indices

x = nan(size(A));
for i = (r+1):(size(A,1)-r)
    x(i, :) = A(i + sample_ixs, :) \ v(i + sample_ixs);
end
x = x';

params_hat = struct();
params_hat.a0 = x(1, :);
params_hat.a1 = sqrt(x(2,:).^2 + x(3,:).^2);
params_hat.p1 = atan2(x(3,:), x(2,:));
params_hat.a2 = sqrt(x(4,:).^2 + x(5,:).^2);
params_hat.p2 = atan2(x(5,:), x(4,:));
params_hat.w1 = w1;
params_hat.w2 = w2;

comps_hat = params2comps(params_hat, t, r=r);

% Reconstruction error: MSE normalized by signal variance, excluding boundaries
ii = (r+1):length(v)-r+1;
err = mean((sum(comps_hat(:,ii),1)-v(ii)').^2)/var(v(ii));
