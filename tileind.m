function n = tileind(varargin)
% Tile index location from row,col
% INPUT:
%	row, col - tile location
%       - negative values start count from last row/col
%       - can be index or one-hot logical vector
%	t - optional, layout object. default = current figure's. 
% SYNTAX:
%   tileind(row,col)
%   tileind(row,col,t)

assert(nargin==2 || nargin==3);

row = varargin{end-1};
col = varargin{end};

if nargin == 3
    t = varargin{1};
elseif nargin == 2
    if isempty(findall(0,'type','figure'))
        error("Cannot find tiled layout - There are no open figures");
    end
    t = get(gcf(),"Children");
    if length(t) ~= 1 || ~contains(class(t),"TiledChartLayout")
        error("Current figure does not have a tiled layout");
    end
end

nrows = t.GridSize(1);
ncols = t.GridSize(2);

if length(row) ~= 1
    assert(length(row) == nrows && sum(row) == 1, ...
        "row location should be either an index, or one-hot logical vector with length==[number of grid rows]");
    row = find(row);
end
if length(col) ~= 1
    assert(length(col) == ncols && sum(col) == 1, ...
        "col location should be either an index, or one-hot logical vector with length==[number of grid cols]");
    col = find(col);
end

if row < 1, row = nrows - row; end
if col < 1, col = ncols - col; end
assert(row <= nrows && col <= ncols, sprintf("Tile location (%d,%d) exceeds layout (%d,%d)", row, col, nrows, ncols));
n = (row - 1) * ncols + col;
