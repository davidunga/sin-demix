function plot_scaled(varargin)

figure();
hold on;
for i = 1 : length(varargin)
    x = varargin{i};
    x = x - min(x);
    x = x / max(x);
    plot(x);
end