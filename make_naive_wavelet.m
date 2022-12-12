function [w,t] = make_naive_wavelet(Fs,f,dur)

t = -.5*dur:(1/Fs):.5*dur;
w = cos(2*pi*f*t);