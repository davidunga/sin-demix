function [comps,opt] = wt_demix2(v,Fs,ws)

df = 4;

has_dc = ws(1)==0;
ws = ws(ws>0);

[WT, WFrqs]=cwt(v,Fs,VoicesPerOctave=12,TimeBandwidth=60);

wght = (ws/(2*pi)).^2;
wght = wght/sum(wght);
wght = wght/(1/length(ws));

comps = nan([length(ws)+1, length(v)]);
for i = 1 : length(ws)
    f = ws(i)/(2*pi);
    ii = WFrqs >= (f-df) & WFrqs <= (f+df);
    ww = WT(ii,:);
    [~,mxi] = max(sum(abs(ww),2));
    comps(i+1,:) = real(ww(mxi,:)).*wght(i);
end

if has_dc
    comps(1,:) = v - sum(comps(2:end,:),1);
else
    comps = comps(2:end,:);
end

opt.r = 100;