function xx = morlet_tform(x,Fs,f,opts)

arguments
    x
    Fs
    f
    opts.n=1
    opts.sigmas=4
end

opts = namedargs2cell(opts);
w = make_morlet(Fs,f,opts{:});


r = floor(length(w)/2)+1;

nConv = length(x) + length(w) - 1;

X = fft(x,nConv);

% create wavelet
W = fft(w,nConv);
W = W./max(abs(W)); % normalize

% convolve
xx = ifft(W.*X);
xx = xx(r:r+length(x)-1);

