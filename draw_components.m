function draw_components(Fs,comps,v)
% Fs - sampling freq
% comps - components
% v - demixed (original) signal, optional.

if nargin==2
    v=[];
end

t = (0:(size(comps,2)-1))'/Fs;

figure();
tl = tiledlayout(4,1);

nexttile(1);
hold on;
if ~isempty(v)
    plot(t,v,'-r',LineWidth=2,DisplayName="Original");
end
plot(t,sum(comps,1),'-b',DisplayName="Reconstructed");
plot(t,sum(comps(1,:),1),'-c',DisplayName="DC");
legend(Location="best",Orientation="horizontal");
title("Signal"); xlim(t([1,end]));

nexttile(2); plot(t,comps(1,:),'c-'); title("DC"); xlim(t([1,end]));
nexttile(3); plot(t,comps(2,:),'-'); title("Sin1"); xlim(t([1,end]));
nexttile(4); plot(t,comps(3,:),'-'); title("Sin2"); xlabel("t"); xlim(t([1,end]));

linkaxes(findobj(tl,'Type','Axes'),'x');
