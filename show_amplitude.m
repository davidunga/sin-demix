function show_amplitude(x,n)

if nargin==1, n=1; end

[a,u,l,s,ui,li] = amplitude(x,n);

opts = {'LineWidth',1};

figure();
hold on;
plot(s,'g-',opts{:});
plot(x,'k-',opts{:});
plot(u,'r-',opts{:});
plot(l,'b-',opts{:});
