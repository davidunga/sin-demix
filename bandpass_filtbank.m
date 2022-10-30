function fb = bandpass_filtbank(N,Fs,f,wv)
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
%       parameters for morse wavelet- [gamma, powerBandWidth]

if ~exist("wv","var")
    wv = "morse";
end

vpo = 48;
params = struct( ...
    SignalLength=N, ...
    SamplingFrequency=Fs, ...
    FrequencyLimits=2.^(log2(f) + (1/vpo) * [-1.0001,1]), ...
    VoicesPerOctave=vpo);

if isnumeric(wv)
    assert(length(wv)==2);
    wv(2) = morseBandWidth_Power2Time(f,wv(2));
    if wv(2)<wv(1) || wv(2)>40*wv(1)
        error("Invalid Params: gamma=%2.2f, TBW=%2.2f",wv(1),wv(2));
    end
    params.WaveletParameters = wv;
else
    params.Wavelet = wv;
end

params = namedargs2cell(params);
fb = cwtfilterbank(params{:});
