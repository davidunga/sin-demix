function show_extremas(x)

figure();
hold on;
plot(x,'g');

locmax = islocalmax(x);
locmin = islocalmax(-x);

if max(length(locmin),length(locmax)) > 30
    marker = ".";
else
    marker = "o";
end

plot(find(locmax),x(locmax),"r" + marker);
plot(find(locmin),x(locmin),"b" + marker);