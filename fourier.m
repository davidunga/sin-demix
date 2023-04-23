function [Fx,frq] = fourier(x, Fs, opts)
% Single-sided Fourier transform

arguments
    x
    Fs
    opts.show = false
end

Fx = fft(x);
n = ceil(length(x)/2) + 1;
frq = (0:(n-1)) / length(x) * Fs;
Fx = Fx(1:n) / length(x);
Fx(2:end) = 2 * Fx(2:end);

if opts.show
    FX = abs(Fx);
    figure();
    tiledlayout(2,1,TileSpacing="compact",Padding="compact");
    nexttile();
    plot((0:(length(x)-1))/Fs,x,'g');
    ylabel('Signal');
    xlabel('t [sec]');
    nexttile();
    plot(frq,FX);
    xlabel('Frq [Hz]');
    ylabel('Amp');
    grid minor;
end