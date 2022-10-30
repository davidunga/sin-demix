function [WTs,fb] = wtbandpass(x,Fs,f,wv)

if ~exist("wv","var"), wv="morse"; end

WTs = nan([length(f),length(x)]);
for i = 1 : length(f)
    fb(i) = bandpass_filtbank(length(x),Fs,f(i),wv);
    WT = cwt(x,FilterBank=fb(i));
    assert(size(WT,1)==3);
    WTs(i,:) = WT(2,:);
end
