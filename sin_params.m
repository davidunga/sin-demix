function [a_new, p_new] = sin_params(x, t,w)
Fs=length(t)/(t(end)-t(1));

T = 2*pi/w;
wfrq = w/(2*pi);
lgw = log2(w/(2*pi));
frqlim = 2.^[floor(lgw), ceil(lgw)];
[wt,frqs] = cwt(x,Fs,FrequencyLimits=frqlim);
[~,row] = max(sum(abs(wt),2));
assert(frqs(row+1) < wfrq);
assert(frqs(row-1) > wfrq);
a_new = abs(wt(row, :));
p_new = angle(wt(row, :));

ifm = round(T * Fs);
ito = length(x) - round(T * Fs);

a_new(1:ifm) = a_new(ifm);
a_new(ito:end) = a_new(ito);

p_new(1:ifm) = p_new(ifm);
p_new(ito:end) = p_new(ito);
return;

mx_ixs = find((x(2:end-1) > x(1:end-2)) & x(2:end-1) > x(3:end)) + 1;
mn_ixs = find((x(2:end-1) < x(1:end-2)) & x(2:end-1) < x(3:end)) + 1;
if mn_ixs(1) < mx_ixs(1)
    mn_ixs = mn_ixs(2:end);
end
assert(mn_ixs(1) > mx_ixs(1));
ixs = sort([mx_ixs, mn_ixs]);
ifm = ixs(1);
ito = ixs(end);

a_new = nan(size(x));
a_new(ifm:ito) = interp1(ixs, smoothdata(abs(x(ixs)), "movmean", 5), ifm:ito,"linear");
a_new(1:ifm) = a_new(ifm);
a_new(ito:end) = a_new(ito);


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
fprintf('T/4=%2.2f\n',T/4);
tt = T/2*(0:(length(ixs)-1)) + T/4;
df = t(ixs)' - tt;
df(:) = df(1);
fprintf('Actual offset =%2.2f\n',t(ixs(1)));
%df = mod(df, T);
%p_new1 = mod(t(ixs),T/2)-T/4;
p_new1 = -df*w;

% close all
% figure();
% plot(t, x, 'k');
% hold on;
% plot(t(ixs),x(ixs),'ro');
% plot(t, a_new(:).*sin(w.*t(:)), 'b');
% 

p_new = nan(size(x));
p_new(ifm:ito) = interp1(ixs, p_new1, ifm:ito,"linear");
p_new(1:ifm) = p_new(ifm);
p_new(ito:end) = p_new(ito);
