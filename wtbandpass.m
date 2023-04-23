function [WTs, fb] = wtbandpass(x,Fs,f,opts)

arguments
    x
    Fs
    f
    opts.wv = "morse"
    opts.interp = true
end

WTs = nan([length(f),length(x)]);
scales = nan([length(f),length(x)]);
for i = 1 : length(f)
    fb(i) = bandpass_filtbank(length(x),Fs,f(i),opts.wv);
    [WT,WT_frqs,~,~,scales(i,:)] = cwt(x,FilterBank=fb(i));
    assert(size(WT,1)==3);

    if opts.interp
    end
    
    WTs(i,:) = WT(2,:);
end