function [comps,opt] = linear_demix_fixedphase(v, Fs, ws, ph)
% Linearly demix signal [v] into a sum of sin-like time varying components with
% angular fequencues [ws], i.e., v is assumed to be a mixture of the form:
%   v(t) = C1(t) + C2(t) + .. Cn(t)
% Where Ci(t) is an amplitude and phase varying sin(), with frequnecy ws(i):
%   Ci(t) = amp(t) * sin(ws(i)*t + phase(t))
% Or, if ws(i) = 0: Ci(t) = amp(t).
% INPUT:    see arguments.
% OUTPUT:   comps - [length(ws),length(v)] array of components, comps(i,:) = Ci(t)
%           opt - struct of option parameters used

arguments
    v                   % 1d signal (vector)
    Fs                  % v's sampling frequency
    ws
    ph
end

assert(all(ws>0));

v = v(:);
t = (0:(length(v)-1))'/Fs;
A = ones([length(t), length(ws)+1]);
for j = 1:length(ws)
    A(:, j+1) = sin(ws(j)*t + ph(j));
end

r_dur = 1*pi/min(ws);
n_samples = 2 * round(r_dur*max(ws)/2*pi) + 1;

r = round(r_dur*Fs);
sample_ixs = round(linspace(-r, r, n_samples));

x = nan(size(A));
for i = (r+1):(size(A,1)-r)
    x(i, :) = A(i + sample_ixs, :) \ v(i + sample_ixs);
end
x(1:r,:) = repmat(x(r+1,:),[r,1]);
x(end-r+1:end,:) = repmat(x(end-r,:),[r,1]);
assert(~any(isnan(x(:))));

x = x';
A = A';

comps = x.*A;
comps(1,:) = v' - sum(comps(2:end,:),1);

opt = struct();
opt.r = r;
opt.a = x(2:end,:);
opt.p = repmat(ph(:), [1,length(v)]);
opt.ph = ph;