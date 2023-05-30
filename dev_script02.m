
% ------------------------

PLOT_ALL_GS = 0;
PLOT_DMX = 1;

% ------------------------

ret = run_script('ModelGeGi_measure_For_paper_Final3_01__2021_fig2and3__JUNK_005');
ret.GTOT = ret.GE+ret.GI+ret.GL;
res = ret.res;

% ------------------------

if PLOT_ALL_GS
    gs_to_plot = ["gtot", "ge", "gi"];
else
    gs_to_plot = ["gtot"];
end

figure(Position=[300,300,800,length(gs_to_plot)*400]);
tl = tiledlayout(length(gs_to_plot)+2*PLOT_DMX,1,TileSpacing="compact",Padding="compact");
for name = gs_to_plot
    nexttile(); hold on;
    plot(ret.T,ret.(upper(name)),"r",DisplayName="True");
    plot(ret.T,res.(name),"k",DisplayName="Predicted");
    xlabel("time");
    ylabel(name);
    grid on;
    add_legend(title=upper(name), Location="NW");
end

if PLOT_DMX
    nexttile(); hold on;
    plot(ret.T,res.Vdmx.a(1,:),"b",DisplayName="a1");
    plot(ret.T,res.Vdmx.a(2,:),"g",DisplayName="a2");
    xlabel("time");
    ylabel("Amp");
    grid on;
    add_legend(title="Estimated V Amps", Location="NW");

    nexttile(); hold on;
    plot(ret.T,res.Vdmx.p(1,:),"b",DisplayName="p1");
    plot(ret.T,res.Vdmx.p(2,:),"g",DisplayName="p2");
    xlabel("time");
    ylabel("Phase");
    grid on;
    add_legend(title="Estimated V Phases", Location="NW");
end

linkaxes(allaxes(),"x");
xlim([5/100, ret.T(end)-5/100]);

fs = round(ret.expr_info.ws/(2*pi));
fig_name = sprintf("%d %d Hz", fs(1), fs(2));
title(tl,fig_name);


% saveas(gcf, fig_name + ".jpg");