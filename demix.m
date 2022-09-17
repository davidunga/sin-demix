function result = demix(v, t, w1, w2, params)
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

v = v(:);
t = t(:);

A = [ones(size(t)), sin(w1*t), cos(w1*t), sin(w2*t), cos(w2*t)];
if exist('params', 'var')
    A(:,2) = A(:,2) .* params.a1(:);
    A(:,3) = A(:,3) .* params.a1(:);
    A(:,4) = A(:,4) .* params.a2(:);
    A(:,5) = A(:,5) .* params.a2(:);
end

x = nan(size(A));
for i = 3:(size(A,1)-2)
    ii = i + (-2:2);
    x(i, :) = A(ii, :) \ v(ii);
end

x(1:2,:) = x([3,3],:);
x(end-1:end,:) = x([end-2, end-2],:);
x = x';

result = struct();
result.a0 = x(1, :);
result.a1 = sin_amp(x(2,:), x(3,:), w1, t);
result.p1 = atan2(x(3,:), x(2,:));
result.a2 = sin_amp(x(4,:), x(5,:), w2, t);
result.p2 = atan2(x(5,:), x(4,:));
result.w1 = w1;
result.w2 = w2;
% %return;
% if exist('params', 'var')
%     result.a1 = params.a1;
%     result.a2 = params.a2;
%     %result.p1 = .5 * (result.p1 + params.p1);
%     %result.p2 = .5 * (result.p2 + params.p2);
% else
%     result.a1 = sin_amp(x(2,:), x(3,:), w1, t);
%     result.a2 = sin_amp(x(4,:), x(5,:), w2, t);
%     result = demix(v, t, w1, w2, result);
% end

function a_new = sin_amp(x1, x2, w, t)
a = sqrt(x1.^2 + x2.^2);
a_new = a;
return;
p = atan2(x2, x1);
x = a .* sin(w.*t' + p);
mx_ixs = find((x(2:end-1) > x(1:end-2)) & x(2:end-1) > x(3:end)) + 1;
mn_ixs = find((x(2:end-1) < x(1:end-2)) & x(2:end-1) < x(3:end)) + 1;
assert(mn_ixs(1) > mx_ixs(1));
ixs = sort([mx_ixs, mn_ixs]);
ifm = ixs(1);
ito = ixs(end);

a_new = nan(size(x));
a_new(ifm:ito) = interp1(ixs, abs(x(ixs)), ifm:ito,"linear");
a_new(1:ifm) = a_new(ifm);
a_new(ito:end) = a_new(ito);
%a_new = smoothdata(a_new, "movmean", 5);

function [a_new, p_new] = reconstruct_sin(x1, x2, w, t)
a = sqrt(x1.^2 + x2.^2);
p = atan2(x2, x1);
x = a .* sin(w.*t' + p);

%[a_new, p_new] = sin_params(x, t,w);

% p_new = nan(size(a_new));
% ii = abs(x1) < abs(x2);
% p_new(ii) = acos(x1(ii)./a_new(ii));
% p_new(~ii) = asin(x2(~ii)./a_new(~ii));
% p_new(:) = median(p_new);
% %ii = abs(a_new) > abs(x1);
% %p_new = acos(x1./a_new);
% return;

mx_ixs = find((x(2:end-1) > x(1:end-2)) & x(2:end-1) > x(3:end)) + 1;
mn_ixs = find((x(2:end-1) < x(1:end-2)) & x(2:end-1) < x(3:end)) + 1;
if mn_ixs(1) < mx_ixs(1)
    mn_ixs = mn_ixs(2:end);
end
assert(mn_ixs(1) > mx_ixs(1));
ixs = sort([mx_ixs, mn_ixs]);
ifm = ixs(1);
ito = ixs(end);


dbg = false;
if dbg
    figure();
    plot(x, 'k'); hold on; plot(ixs,x(ixs),'ro');
    figure();
    plot(abs(x(ixs)));
end

a_new = nan(size(x));
a_new(ifm:ito) = interp1(ixs, smoothdata(abs(x(ixs)), "movmean", 3), ifm:ito,"linear");
a_new(1:ifm) = a_new(ifm);
a_new(ito:end) = a_new(ito);

T = 2*pi/w;
p0 = t(ixs(1)) - T/4;
df = smoothdata(diff(t(ixs)'),"movmean", 5) - T/2;
p_new1 = [p0, p0 + cumsum(df)];

p_new = nan(size(x));
p_new(ifm:ito) = interp1(ixs, p_new1, ifm:ito,"linear");
p_new(1:ifm) = p_new(ifm);
p_new(ito:end) = p_new(ito);
p_new = -p_new*w;
return;
% Fs = 120;
% t = 0:(1/Fs):20;
% w=3;
% x = sin(w*t + 000000000.54);
% mx_ixs = find((x(2:end-1) > x(1:end-2)) & x(2:end-1) > x(3:end)) + 1;
% 
% tt = t(mx_ixs);
% T = 2*pi/w;
% m = mod(tt - T/4, T);
% m(m > .99*T) = 0;
% 

T = 2*pi/w;

%df = t(mx_ixs)' - T*(0:(length(mx_ixs)-1)) - T/4;
%p_new1 = df;

tt = T/2*(0:(length(ixs)-1)) + T/4;
df = t(ixs)' - tt;
%df = mod(df, T);
%p_new1 = mod(t(ixs),T/2)-T/4;
p_new1 = -df*w;


dbg = false;
if dbg
    figure();
    plot(t, x, 'k'); hold on; plot(t(ixs),x(ixs),'ro');
    for i = 1 : length(tt)
        [~,mni] = min(abs(tt-t(i)));
        plot(tt(i),x(mni),'m*');
    end
end

p_new = nan(size(x));
p_new(ifm:ito) = interp1(ixs, p_new1, ifm:ito,"linear");
p_new(1:ifm) = p_new(ifm);
p_new(ito:end) = p_new(ito);



