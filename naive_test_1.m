
Fs = 5000;
dur=3;
t = (0:(dur*Fs-1))/Fs;

f1 = 15;
f2 = 1.9*f1;

L = length(t);
ix1 = round(.2*L):round(.6*L);
ix2 = round(.4*L):round(.8*L);
%ix2 = round(.5*L):round(.8*L);

C1 = zeros(size(t));
C1(ix1) = sin(2*pi*f1*t(ix1));

C2 = zeros(size(t));
C2(ix2) = sin(2*pi*f2*t(ix2));

gt_comps = [C1;C2];
v = sum(gt_comps,1);

T = min(f1,f2)^-1;
[wv1,t1] = make_naive_wavelet(Fs,f1,T);
[wv2,t2] = make_naive_wavelet(Fs,f2,T);

% Check magntiudes:
% Should be = integral{sin(wt)^2}dt t=[-T/2,T/2] = T/2, upto Fs
mag1 = max(conv(sin(2*pi*f1*t),wv1,"same"));
mag2 = max(conv(sin(2*pi*f2*t),wv2,"same"));
fprintf("Wavelet1 magnitude: Analytic=%2.2f Computed=%2.2f\n",.5*(T*Fs+1),mag1);
fprintf("Wavelet2 magnitude: Analytic=%2.2f Computed=%2.2f\n",.5*(T*Fs+1),mag2);

% normalize wavelets:
wv1 = wv1 / mag1;
wv2 = wv2 / mag2;

c1 = conv(v,wv1,"same");
c2 = conv(v,wv2,"same");
comps = [c1;c2];

[comps2,~] = naive_demix(v,Fs,2*pi*[0,f1,f2]);
comps = comps2(2:3,:);

figure();
tiledlayout(size(comps,1),1);
for i=1:size(comps,1)
    nexttile();
    plot(t,gt_comps(i,:),'r'); hold on; plot(t,comps(i,:),'b');
end
legend({'gt','est'});
