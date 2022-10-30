function x = set_amplitude(x,a,n)

if nargin<3
    n=1;
end

[a_curr,~,~,s] = amplitude(x,n);
x = (x - s) ./ a_curr .* a + s;