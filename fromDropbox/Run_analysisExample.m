

%Hi David
%run this script to get the analysis to get ge and gi.
%to change from my analysis to yours change this global  variable:'DavidData'
load VI_example
global DavidData Stable_demix_op;
DavidData = 1; % 0 or 1. If 0 the bandpass approach  and if 1 David's approach

Stable_demix_op = 1;
vl = -0.07;
revs = [0.00 vl -0.08]; % [ve vl vi]
FreqArray= [1 2];
FiltType = 1; %1 for bandpass, 0 for filtfilt 3 for wavelet
hybridCe = 0;
BoostCe = 0;
removeHighLow = 0; %aug 2022 new way to remove artifact
FILTP = [30 0.0001 3 0.83];
compensLowRin = 0; %to prevent z2>z1, lowering the real z2
c = 0;

openFigs = allchild(0);
[ge,gi,gl,re,VC,GT,Zt,cmm,Xfound,ff,ff2,g1,g2,z1,z2,vl,Iin,Iex]=  findGeGi_MultiFreq_v005_temp(V,Iinj,1/dt,c,revs,[0.4 0.7],1,FILTP,-1,BoostCe,0,FreqArray,FiltType,hybridCe,removeHighLow,compensLowRin);% dsds ?????%


close(setdiff(allchild(0),openFigs));
figure();
tiledlayout(3,1);
T = (1:length(ge))*dt;
tlim = [.7,1.3];
ii = T>tlim(1) & T<tlim(2);
nexttile();
plot(T(ii),GE(ii)); hold on; plot(T(ii),ge(ii)); ylabel('GE');
legend({'real','measured'},'Location','best');
nexttile();
plot(T(ii),GI(ii)); hold on; plot(T(ii),gi(ii)); ylabel('GI');
nexttile();
hold on; plot(real(z1)); plot(real(z2));
xlim([2.3417, 3.8995] * 1e4); ylim([2.8787, 3.5353] * 1e7); ylabel('Z');


