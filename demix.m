function result = demix(v, t, w1, w2, smooth)
% Demix time sequence into a time-modulated sum of sins.
% Specifically, for a signal v(t), assumed to be a mixture of the form:
%   v(t) = a0(t) + a1(t)*sin(w1*t + p1(t)) + a2(t)*sin(w2*t + p2(t))
% Finds the time-varying amplitudes {a0,a1,a2} and phases {p1,p2}, under
% the assumption that they are varying slower than max(w1,w2).
% INPUT:
%   v, t - signal and time vectors to demix
%   w1, w2 - frequency of the underlying sins
%   smooth - bool, smooth params between pointwise solutions? default=true
% OUTPUT:
%   struct with esimated parameters.

if ~exist('smooth', 'var'), smooth = true; end

v = v(:);
t = t(:);

A = [ones(size(t)), sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t)];
x = nan(size(A));
for i = 3:(size(A,1)-2)
    ii = i + (-2:2);
    x(i, :) = A(ii, :) \ v(ii);
end

x(1:2,:) = x([3,3],:);
x(end-1:end,:) = x([end-2, end-2],:);
x = x';

if smooth
    x = smoothdata(x, 2, "movmean", 2 * 5 + 1);
end

result = struct();
result.a0 = x(1, :);
result.a1 = sqrt(x(2,:).^2 + x(3,:).^2);
result.p1 = atan2(x(3,:), x(2,:));
result.a2 = sqrt(x(4,:).^2 + x(5,:).^2);
result.p2 = atan2(x(5,:), x(4,:));
result.w1 = w1;
result.w2 = w2;
