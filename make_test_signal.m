function [s,t] = make_test_signal(Fs, dur, offset)

gauss_risetime = .2;

% ----
sd = gauss_risetime / 5;
step_dur = gauss_risetime * .75;
delay = gauss_risetime;
% ----

t = fs2tm(Fs, dur);
s = zeros(size(t));

ifm = round((delay + offset) * Fs);
ito = ifm + round(step_dur * Fs);
s(ifm:ito) = 1;

tt = max(t(s>1e-10));
mu = tt + step_dur + gauss_risetime;
g = exp(-((t-mu)/(2*sd)).^2);
s = s + g / max(g);
