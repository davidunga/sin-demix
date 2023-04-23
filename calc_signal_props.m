function [ws_refined, ph, rest_ixs, dmx] = calc_signal_props(x, Fs, ws)

if nargin==2
    ws = sort(dominant_freqs(x,Fs))*2*pi;
end

dmx = wt_demix(x,Fs,ws);
rest_ixs = get_low_variance_window(dmx.dc,Fs,ws);
t = dmx.t;
p = dmx.p + pi/2;

ws_refined = zeros(size(ws));
ph = zeros(size(ws));
for i = 1:length(ws)
    
    % naive estimate = median of d(phase)/dt:
    w_est = median(diff(p(i,rest_ixs)))*Fs;

    % RANSAC, using MAD from naive estimate as margin:
    phase_diff = abs(p(i,rest_ixs) - ws(i)*t(rest_ixs));
    phase_diff_MAD = median(abs(phase_diff-median(phase_diff)));
    ransac_result = ransac_line([t(rest_ixs)',p(i,rest_ixs)'], phase_diff_MAD);

    ws_refined(i) = ransac_result(1);   % omega
    ph(i) = ransac_result(2);   % base angular offset

end