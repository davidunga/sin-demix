
% ------------------------

PLOT_ALL_GS = 0;

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
tl = tiledlayout(length(gs_to_plot),1,TileSpacing="compact",Padding="compact");
for name = gs_to_plot
    disp(name);
    nexttile(); hold on;
    plot(ret.T,ret.(upper(name)),"r",DisplayName="True");
    plot(ret.T,res.(name),"b",DisplayName="Predicted");
    xlabel("time");
    ylabel(name);
    grid on;
    lgnd=legend();
    lgnd.Orientation="horizontal";
    lgnd.Location = "north";
    lgnd.Title.Visible = true;
    lgnd.Title.String = upper(name);
end

fs=round(ret.expr.ws/(2*pi));
fig_name = sprintf("%d %d Hz", fs(1), fs(2));
title(tl,fig_name);
saveas(gcf, fig_name + ".jpg");