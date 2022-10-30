function [WTs,fb] = wtbandpass(x,Fs,f,wv,normP2)

if ~exist("wv","var"), wv="morse"; end
if ~exist("normP2","var")
    normP2=false;
elseif normP2
    assert(isnumeric(wv));
end

fref = {};
if isnumeric(wv) && normP2
    fref = {max(f)};
end

WTs = nan([length(f),length(x)]);
for i = 1 : length(f)
    fb(i) = bandpass_filtbank(length(x),Fs,f(i),wv,fref{:});
    WT = cwt(x,FilterBank=fb(i));
    assert(size(WT,1)==3);
    WTs(i,:) = WT(2,:);
end
