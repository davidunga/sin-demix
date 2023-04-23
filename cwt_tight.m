function ret = cwt_tight(x,Fs,frqs,opts)

arguments
    x
    Fs
    frqs
    opts.wv = "morse"
    opts.r = 1
end

ret = struct();
for i = 1:length(frqs)
    ret(i).fb = bandpass_filtbank(length(x), Fs, frqs(i), wv=opts.wv, r=opts.r);
    [ret(i).WT, ret(i).frqs] = cwt(x,FilterBank=ret(i).fb);
    assert(size(ret(i).WT,1)==2*opts.r+1);
end