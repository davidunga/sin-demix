

Fs = 5000;
dur = 3;
t = (0:(dur*Fs-1))/Fs;

f1 = 51;
f2 = 82;

L = length(t);
ix1 = round(.2*L):round(.6*L);
ix2 = round(.5*L):round(.8*L);

C1 = zeros(size(t));
C1(ix1) = sin(2*pi*f1*t(ix1));

C2 = zeros(size(t));
C2(ix2) = sin(2*pi*f2*t(ix2));

x = C1+C2;

n=4;
sigmas=5;
c1 = morlet_tform(x,Fs,f1,n=n,sigmas=sigmas);
c2 = morlet_tform(x,Fs,f2,n=n,sigmas=sigmas);

figure();
tiledlayout(3,1,TileSpacing="compact",Padding="compact");
nexttile();
plot(t,x,'r',t,real(c1)+real(c2),'b');
nexttile();
plot(t,C1,'r',t,real(c1),'b');
nexttile();
plot(t,C2,'r',t,real(c2),'b');


