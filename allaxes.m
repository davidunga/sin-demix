function axs = allaxes()

try
    h = gctl();
catch
    h = gcf();
end
axs = findobj(get(h,'Children'),'Type','Axes');