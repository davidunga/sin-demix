function tbw = morseBandWidth_Power2Time(f,halfpbw)
% morse wavelet time bandwidth parameter from half power bandwidth around
% frequency

% reference values (arbitrarily chosen):
% for 100Hz, and time bandwidth 10, the half power bandwidth is 52.34:
REF_F = 100;
REF_TBW = 10;
REF_HPBW = 52.34;

factor = (REF_TBW * REF_HPBW^2) * (f/REF_F)^2;

tbw = factor * halfpbw^-2;