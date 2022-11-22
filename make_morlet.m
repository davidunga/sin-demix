function w = make_morlet(Fs,f,opts)

arguments
    Fs
    f
    opts.n = 4
    opts.sigmas = 4
end

sig = opts.n/(2*pi*f);
t = -opts.sigmas*sig:(1/Fs):opts.sigmas*sig;
w = exp(-.5 * t.^2 / sig^2) .* exp(1i*2*pi*f*t);
w = 2 * w ./ sum(abs(w));

DBG_PLOT = 0;
if DBG_PLOT
    figure();
    hold on;
    plot(t,real(w),'b');
    plot(t,imag(w),'r');
    plot(t,abs(w),'k');
end