
ModelGeGi_measure_For_paper_Final3_01__2021_fig2and3__JUNK_005;

%feb 2023
close all

figure; plot(GE+GI+GL); hold on; plot(ge+gi+gl); legend(["GT","Pred"]); title("Sum of G's");    
hg=gca;
figure; plot(real(z1)); hold on; plot(real(z2)); legend(["z1", "z2"]); title("Z's");
hz=gca;

gsum = gi+ge+gl;
ii=islocalmax(gsum,"MinSeparation",10,"MaxNumExtrema",3);
axes(hg);
plot(find(ii),gsum(ii),'o');

axes(hz);
xline(find(ii));