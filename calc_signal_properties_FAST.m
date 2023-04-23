function ret = calc_signal_properties_FAST(v, Fs)
% Calc omegas and margin of influence, using fft
domfrqs = dominant_freqs(v,Fs);
ret.ws = sort(domfrqs)*2*pi;
ret.margin = round(2*2*pi/min(ret.ws)*Fs);  % margin ~ two periods
return