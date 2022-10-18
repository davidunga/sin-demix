function run_example(seed)

if ~exist('seed', 'var'), seed=randi(1000); end

rng(seed);

% --
% build signal:

w1 = 151 * 2 * pi;
w2 = 282 * 2 * pi;

Fs = 5e3; % sampling freq [Hz]
T = .25; % signal duration [sec]
t = (0:(T*Fs-1))/Fs;

assert(Fs > 2 * max(w1, w2));

% build parameters as baseline values which are (slowly) modulated using
% random sin, defined by random frequency and relative amplitute from these ranges: 
w_rng = min([w1, w2]) / 10 * [.5, 1];
amp_rng = [.1, .2];
phase_amp = 5 * pi / 180;

params = struct();
params.a0 = sin_modulate(t, 0.1, randrng(w_rng), randrng(amp_rng));
params.a1 = sin_modulate(t, 0.3, randrng(w_rng), randrng(amp_rng));
params.a2 = sin_modulate(t, 0.2, randrng(w_rng), randrng(amp_rng));
params.p1 = pi/8 + phase_amp * sin(randrng(w_rng) * t);
params.p2 = pi/4 + phase_amp * sin(randrng(w_rng) * t);
params.w1 = w1;
params.w2 = w2;

comps = params2comps(params, t);
v = sum(comps,1); % the mixed signal

% --
% estimate components from mixed signal

comps_hat = demix(v, Fs, w1, w2);

% --
% show results

disp("randseed " + num2str(seed));
evaluate_prediction(comps, comps_hat);
plot_gt_and_prediction(t, comps, comps_hat);

