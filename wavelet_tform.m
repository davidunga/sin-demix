function xx = wavelet_tform(x,Fs,f,opts)

arguments
    x
    Fs
    f
    opts.n=1
    opts.sigmas=4
end

r = floor(length(wv)/2)+1;

nConv = length(x) + length(wv) - 1;

X = fft(x,nConv);

% create wavelet
W = fft(w,nConv);
W = W./max(abs(W)); % normalize

% convolve
xx = ifft(W.*X);
xx = xx(r:r+length(x)-1);

