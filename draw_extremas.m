function draw_extremas(x,t)

if nargin==1
    t = 1:length(x);
end

locmax = find(islocalmax(x));
locmin = find(islocalmax(-x));

figure();
hold on;
plot(t,x,'g');
plot(t(locmax),x(locmax),'r.');
plot(t(locmin),x(locmin),'b.');