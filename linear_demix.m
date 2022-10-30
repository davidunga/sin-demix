function [comps,opt] = linear_demix(v, Fs, ws, options)
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
    ws                  % vector of frequencies, in accending order
    options.eqs = 0     % number of equations to use, default = number of variables
    options.r_dur = 0   % radius of sampling window, in seconds. default = shortest period
end

assert(length(ws)==1 || all(diff(ws) > 0), "Frequencies must be monotonically increasing");

v = v(:);
t = (0:(length(v)-1))'/Fs;
num_dcs = double(ws(1) == 0);
num_vars = 2 * length(ws) - num_dcs;

A = ones([length(t), num_vars]);
for k = 1:(length(ws) - num_dcs)
    w = ws(num_dcs + k);
    i = num_dcs + 2 * k - 1;
    A(:,i:i+1) = [sin(w*t), cos(w*t)];
end

if options.eqs > 0
    eqs = options.eqs;
    assert(eqs >= num_vars, sprintf("Need at least %d equations, got: %d", num_vars, eqs));
else
    eqs = num_vars;
end

if options.r_dur > 0
    r_dur = options.r_dur;
else
    r_dur = 2*pi/max(ws);
end

r = round(r_dur * Fs);  % sampling radius
sample_ixs = round(linspace(-r,r,eqs)); % 0-centered sampling indices

x = nan(size(A));
for i = (r+1):(size(A,1)-r)
    x(i, :) = A(i + sample_ixs, :) \ v(i + sample_ixs);
end
x(1:r,:) = repmat(x(r+1,:),[r,1]);
x(end-r+1:end,:) = repmat(x(end-r,:),[r,1]);

comps = nan([length(t),length(ws) - num_dcs]);
if num_dcs==1
    comps(:,1) = x(:,1);
end
for k = 1:(length(ws)-num_dcs)
    w = ws(num_dcs + k);
    i = num_dcs + 2 * k - 1;
    a = sqrt(sum(x(:,i:i+1).^2,2));
    p = atan2(x(:,i+1), x(:,i));
    comps(:,k + num_dcs) = a .* sin(w .* t + p);
end
comps = comps';

opt.r_dur = r_dur;
opt.r = r;
opt.sample_ixs = sample_ixs;
