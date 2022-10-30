function fb = bandpass_filtbank(N,Fs,f,wv,fref)
% Frequency-localized filter bank.
% Creates a filter bank which includes excactly 3 frequencies: [f-df,f,f+df],
%   where df is the smallest frequency step (="voice") around f.
%
% INPUT:
% N - signal length
% Fs - sampling rate (Hz)
% f - frequency to localize around (Hz)
% wv - wavelet argument.
%   Either:
%       wavelent name- {"morse"} / "amor" / "bump"
%   or:
%       parameters for morse wavelet- [gamma, P2], Where: 1 <= gamma <= P2, P2 <= 40*gamma
%       If P2=-1, P2 is taken to be the largest possible value, depending also on fref.
% fref - reference frequency. can be supplied only if wv is morse parameters.
%   indicates that P2 is in reference to frequency=fref, and should be scaled to fit frequency=f.

if ~exist("wv","var")
    wv = "morse";
end

if exist("fref","var")
    assert(isnumeric(wv));
else
    fref = f;
end

vpo = 48;
params = struct( ...
    SignalLength=N, ...
    SamplingFrequency=Fs, ...
    FrequencyLimits=2.^(log2(f) + (1/vpo) * [-1.0001,1]), ...
    VoicesPerOctave=vpo);

if isnumeric(wv)
    assert(length(wv)==2);
    gamma = wv(1);
    P2 = wv(2);
    if P2 == -1
        P2 = gamma*40;
    end
    params.WaveletParameters = [gamma,P2*(f/fref)^2];
else
    params.Wavelet = wv;
end

params = namedargs2cell(params);
fb = cwtfilterbank(params{:});
