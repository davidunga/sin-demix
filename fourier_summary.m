function fourier_summary(x,Fs,opts)
% Summary of dominant Fourier coeffs

arguments
    x
    Fs
    opts.N = 5
    opts.rel = .05
end

[Fx,frqs] = fourier(x,Fs);
FX = abs(Fx);

ixs = find(islocalmax(FX));

[~,si] = sort(-FX(ixs));
N = min(length(ixs), opts.N);
ixs = ixs(si(1:N));

ixs = ixs(FX(ixs) > opts.rel * FX(ixs(1)));
ixs = [1,ixs];

if length(ixs) == 1
    FXmax = FX(ixs(1));
else
    FXmax = FX(ixs(2));
end

for i = ixs(:)'
    fprintf(" f=%-6.2f  w=%-6.2f  amp=%-6.2f (%2.2f of max)\n", frqs(i), frqs(i) *2 * pi, FX(i), FX(i) / FXmax);
end


