
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
w_rng = min([w1, w2]) / 10 * [.5, 1];
amp_rng = [.1, .2];
max_phase1 = 20*pi/w1;
max_phase2 = 20*pi/w2;

params = struct();
params.a0 = sin_modulate(t, 0.1, randrng(w_rng), randrng(amp_rng));
params.a1 = sin_modulate(t, 0.3, randrng(w_rng), randrng(amp_rng));
params.a2 = sin_modulate(t, 0.2, randrng(w_rng), randrng(amp_rng));
params.p1 = .5 * max_phase1 + .1 * max_phase1 * sin(randrng(w_rng) * t);
params.p2 = .3 * max_phase1 + .1 * max_phase1 * sin(randrng(w_rng) * t);
params.w1 = w1;
params.w2 = w2;

params.p1 = zeros(size(params.p1))+.95*max_phase1;
params.p2 = zeros(size(params.p2))+.95*max_phase2;

fprintf("gt p1 avg= %2.3f\n", mean(params.p1));
fprintf("gt p2 avg= %2.3f\n", mean(params.p2));

%params.p1 = mod(params.p1, 2*pi/params.w1);
%params.p2 = mod(params.p2, 2*pi/params.w2);

[v, comps] = params2signal(params, t); % the mixed signal

% --
% estimate components from mixed signal

params_hat = demix(v, t, w1, w2);
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
