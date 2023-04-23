function plt(varargin)
% plt y1 y2 ...
% plt x=t y1 y2 ...
% plt ax=gca ...

params.x = [];
params.ax = [];
params.xl = [];
params.yl = [];
params.norm = "none";

args_mask = true(size(varargin));
for i = 1 : nargin
    
    if strcmp(varargin{i}, '-h')
        varargin{i} = 'ax=gca';
    end

    tokens = strsplit(varargin{i},'=');
    if length(tokens) == 1
        continue;
    end

    assert(length(tokens)==2);

    key = tokens{1};
    assert(isfield(params,key));
    if strcmp(key, "norm")
        params.(key) = tokens{2};
    else
        params.(key) = evalin("caller", tokens{2});
    end
    
    args_mask(i) = false;
end

varargin = varargin(args_mask);

if isempty(params.ax)
    figure();
    params.ax = gca();
end

hold(params.ax, 'on');

lgnd = {};
for i = 1 : length(varargin)
    y = evalin("caller", varargin{i});
    if ismatrix(y) && size(y,1) > 1 && size(y,2) > 1
        for row = 1 : size(y,1)
            lgnd{end+1} = sprintf("%s row %d",varargin{i},row);
            do_plot(y(row,:),params);
        end
    else
        lgnd{end+1} = varargin{i};
        do_plot(y,params);
    end
end
lgnd=legend(lgnd{:});
lgnd.Interpreter = 'none';
lgnd.Location = 'best';
lgnd.AutoUpdate = 'off';

if ~isempty(params.xl)
    xline(params.xl);
end
if ~isempty(params.yl)
    yline(params.yl);
end

function do_plot(y,params)
if ~isempty(params.x)
    xx = params.x;
else
    xx = 1:length(y);
end
switch params.norm
    case "max"
        y = (y - min(y)) / (max(y) - min(y));
    case "std"
        y = (y - mean(y)) / std(y);
    otherwise
        assert(strcmp(params.norm,"none"));
end
plot(params.ax, xx, y);

