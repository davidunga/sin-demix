global configs
cfg = struct();
i=0;
for fnc = ["linear_demix2", "wt_demix2d"]
    for stable = [false, true]
        for dc_from_v = [false, true]
            for Iclean = [false, true]
                for V0clean = [false, true]

                        i=i+1;
                        cfg(i).ix = i+64;
                        cfg(i).fnc = fnc;
                        cfg(i).stable = stable;
                        cfg(i).dc_from_v = dc_from_v;

                        cfg(i).Iclean = Iclean;
                        cfg(i).V0clean = V0clean;

                        
                        if cfg(i).ix ~= 78
                            continue;
                        end

                        configs = cfg(i);
                        dev_script02;
                        set(gcf,'position',[100,100,1200,600]);
                        name = sprintf("%d %s stable%d dcFromV%d Iclean%d V0clean%d gtotCalc%s", ...
                            i+64, fnc, stable, dc_from_v, Iclean, V0clean, "VAR2");

                        title(name, Interpreter="none");
                        saveas(gcf,name + ".png");

                end
            end
        end
    end
end