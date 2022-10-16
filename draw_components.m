function draw_components(t,comps,v)

if nargin==2
    v=[];
end

figure();
tl = tiledlayout(4,1);

nexttile(); 
hold on;
if ~isempty(v)
    plot(t,v,'-m',LineWidth=2,DisplayName="Original");
end
plot(t,sum(comps,1),'-b',DisplayName="Reconstructed");
plot(t,sum(comps(1,:),1),'-c',DisplayName="DC");
title("Signal");

nexttile(); plot(t,comps(1,:),'-'); title("DC");
nexttile(); plot(t,comps(2,:),'-'); title("Sin1");
nexttile(); plot(t,comps(3,:),'-'); title("Sin2");

linkaxes(allchild(tl),'x');
