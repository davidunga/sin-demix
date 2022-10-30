function compare_components(Fs,compss,optss)

if nargin==2
    names = string(1:length(compss));
else
    names = string(round([optss.r_dur]*1000,1)) + "ms";
end

figure();
tl = tiledlayout(4,1);

h = [];
for i = 1 : length(compss)
    comps = [sum(compss{i},1);compss{i}];
    t = (0:(size(comps,2)-1))'/Fs;
    for c = 1 : 4
        nexttile(c); hold on;
        h(end+1)=plot(t, comps(c,:),DisplayName=names(i));
    end
end
legend(h(1:4:end));

linkaxes(findobj(tl,'Type','Axes'),'x');