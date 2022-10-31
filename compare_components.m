function compare_components(Fs,compss,names)
% compare_components(Fs,{comps1,comps2,..})
% compare_components(Fs,{comps1,comps2,..},["name1","name2",..])

if nargin==2
    names = string(1:length(names));
else
    names = string(names);
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
