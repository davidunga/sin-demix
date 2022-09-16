
clear all;

% --
% build signal:

w1 = 151 * 2 * pi;
w2 = 282 * 2 * pi;

Fs = 5e3; % sampling freq [Hz]
T = .5; % signal duration [sec]
t = (0:(T*Fs-1))/Fs;

% build parameters as baseline values which are (slowly) modulated using
% random sin, defined by random frequency and relative amplitute from these ranges: 
modulation_w = min([w1, w2]) / 10 * [.5, 1];
modulation_amp = [.1, .2];

params = struct();
params.a0 = sin_modulate(0.1, t, modulation_w, modulation_amp);
%params.a0 = piecewise_const(length(t), [-.1, -.05]);
params.a1 = sin_modulate(0.03, t, modulation_w, modulation_amp);
params.a2 = sin_modulate(0.02, t, modulation_w, modulation_amp);
params.p1 = sin_modulate(0.05, t, modulation_w, modulation_amp);
params.p2 = sin_modulate(0.05, t, modulation_w, modulation_amp);
params.w1 = w1;
params.w2 = w2;

[v, comps] = params2signal(params, t); % the mixed signal

% --
% estimate components from mixed signal

params_hat = demix(v, t, w1, w2, false);
[v_hat, comps_hat] = params2signal(params_hat, t);

% --
% show:

figure();
plot_sines(comps, t, true);
plot_sines(comps_hat, t, false);
legend({'GT v', 'GT dc', 'GT sin1 ', 'GT sin2', 'Estimated'});
title('Signal & Components');

figure();
plot_params(params, t, true);
plot_params(params_hat, t, false, true);
title('Parameters');
