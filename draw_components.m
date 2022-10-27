function draw_components(varargin)
% draw_components(comps,Fs[,v])
%   - comps = (3,:) components matrix,
%   - Fs = sampling rate
%   - v (optional) = the original (demixed) signal
%
% draw_components(result[,v])
%   - result = result struct
%   - v (optional) = the original (demixed) signal

v = [];
if isstruct(varargin{1})
    Fs = varargin{1}.Fs;
    comps = [varargin{1}.dc;varargin{1}.sin1;varargin{1}.sin2];
    if length(varargin)==2
        v=varargin{2};
    end
else
    comps = varargin{1};
    Fs = varargin{2};
    if length(varargin)==3
        v=varargin{3};
    end
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
plot(t,comps(1,:),'-c',DisplayName="DC");
legend(Location="best",Orientation="horizontal");
title("Signal"); xlim(t([1,end]));

nexttile(2); plot(t,comps(1,:),'c-'); title("DC"); xlim(t([1,end]));
nexttile(3); plot(t,comps(2,:),'-'); title("Sin1"); xlim(t([1,end]));
nexttile(4); plot(t,comps(3,:),'-'); title("Sin2"); xlabel("t"); xlim(t([1,end]));

linkaxes(findobj(tl,'Type','Axes'),'x');
