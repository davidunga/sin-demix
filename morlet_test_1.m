

Fs = 1000;
dur = 5;
t = (0:(dur*Fs-1))/Fs;

f1 = 151;
f2 = 282;

L = length(t);
ix1 = round(.2*L):round(.6*L);
ix2 = round(.5*L):round(.8*L);

C1 = zeros(size(t));
C1(ix1) = sin(2*pi*f1*t(ix1));

C2 = zeros(size(t));
C2(ix2) = sin(2*pi*f2*t(ix2));

x = C1+C2;

c1 = conv(x,make_morlet(Fs,f1,n=1),'same');
c2 = conv(x,make_morlet(Fs,f2,n=1),'same');


