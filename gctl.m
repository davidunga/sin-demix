function t = gctl()
% Get current tiled layout

if isempty(findall(0,'type','figure'))
    error("There are no open figures");
end
t = get(gcf(),"Children");
if length(t) ~= 1 || ~contains(class(t),"TiledChartLayout")
    error("Current figure does not have a tiled layout");
end