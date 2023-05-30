function lgnd = add_legend(opts)

arguments
    opts.Orientation = "horizontal";
    opts.Location = "best";
    opts.AutoUpdate = false;
    opts.title = []
end

if strlength(opts.Orientation) == 1
    opts.Orientation = replace(lower(opts.Orientation),["h","V"], ["horizontal", "vertical"]);
end

if strlength(opts.Location) <= 3
    opts.Location = replace(lower(opts.Location),["n", "e", "w", "s", "o"], ...
        ["north", "east", "west", "south", "outside"]);
end

lgnd=legend();

lgnd.AutoUpdate = opts.AutoUpdate;
lgnd.Orientation = opts.Orientation;
lgnd.Location = opts.Location;

if ~isempty(opts.title)
    lgnd.Title.Visible = true;
    lgnd.Title.String = opts.title;
end
