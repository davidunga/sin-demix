function [dom_frqs, dom_amps, FT, frqs] = dominant_freqs(x, Fs, opts)
% Find dominant frequencies (ignoring dc component)

arguments
    x
    Fs
    opts.n = Inf            % max number of frequencies
    opts.mindist = 1        % minimal distance between frequency peaks [Hz]
    opts.band = [0, Inf]    % frequency band
    opts.show = false       % show result
    opts.min_prom = .1      % min peak prominence relative to amplitude range
end

% -----------------------------------------------------------------
% Get fourier peaks

[FT, frqs] = fourier(x - mean(x), Fs);
amps = abs(FT);
if opts.mindist > 0
    amps = movmean(amps,opts.mindist/2,SamplePoints=frqs);
end
band_ixs = (frqs >= opts.band(1)) & (frqs <= opts.band(2));
amp_range = max(amps(band_ixs)) - min(amps(band_ixs));
[peak_amps, peak_frqs] = findpeaks(amps(band_ixs), frqs(band_ixs), ...
    MinPeakDistance=opts.mindist, MinPeakProminence=opts.min_prom*amp_range);

if length(peak_amps) > opts.n
    [~,si] = sort(peak_amps, 'descend');
    si = sort(si(1:opts.n));
    peak_amps = peak_amps(si);
    peak_frqs = peak_frqs(si);
end

% -----------------------------------------------------------------
% Refine- replace each peak with nearest spline local maxima

refined_peaks = struct();
spln = spline(frqs, amps);
dspln = fnder(spln);
for i = 1:length(peak_amps)
    ix = find(frqs==peak_frqs(i));
    refined_peaks(i).I = frqs([ix-1, ix+1]);
    fnzero_result = fnzeros(dspln, refined_peaks(i).I);
    refined_peaks(i).frq = fnzero_result(1);
    refined_peaks(i).amp = ppval(spln, refined_peaks(i).frq);
end

dom_amps = [refined_peaks.amp];
dom_frqs = [refined_peaks.frq];

% -----------------------------------------------------------------

if opts.show
    figure();
    tiledlayout(1,1,TileSpacing="compact",Padding="compact");
    nexttile(); hold on;
    plot(frqs,abs(FT),'k'); 
    for i = 1:length(refined_peaks)
        ff = linspace(refined_peaks(i).I(1), refined_peaks(i).I(end), 500);
        plot(ff, ppval(spln, ff), 'c');
        plot(dom_frqs(i), dom_amps(i), 'r*');
    end
    xline([min(frqs(band_ixs)), max(frqs(band_ixs))]);
    xlim([frqs(1), frqs(end)]);
    xlabel('Frq [Hz]');
    ylabel('|FFT|');
    grid minor;
end
