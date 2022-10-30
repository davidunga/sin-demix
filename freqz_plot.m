function freqz_plot(fb,opt)
% wrapper for freqz() plot

arguments
    fb
    opt.ax = []
end

if mean(fb.FrequencyLimits) < 1000
    unitsFactor = 1 / 1000;
    tickstep = 10;
    frqUnits = 'Hz';
else
    unitsFactor = 1;
    tickstep = 100;
    frqUnits = 'KHz';
end

flims = round([fb.FrequencyLimits(1) - 100, fb.FrequencyLimits(2) + 100] / tickstep) * tickstep;

if isempty(opt.ax)
    figure();
else
    axes(opt.ax);
end

freqz(fb);

xtcks = flims(1):tickstep:flims(2);
xticks(xtcks * unitsFactor);
xticklabels(xtcks);
xlim(flims * unitsFactor);
xlabel('Freq [Hz]');

if isempty(opt.ax)
    title(sprintf("[%d, %d]",fb.WaveletParameters));
end
