function [comps, opt] = chatGPT_estimate_parameters(v, Fs, w1, w2)

% Define the model equation
model = @(t, x) x(1) + x(2)*sin(w1*t + x(4)) + x(3)*sin(w2*t + x(5));

% Define the cost function as the sum of squared residuals
cost = @(x) sum((v - model(1:length(v), x)).^2);

% Set initial parameter values
x0 = [mean(v), std(v)/2, std(v)/2, 0, 0];

% Perform optimization to find best parameter values
x = fminsearch(cost, x0);

% Extract parameter values
d = x(1)*ones(size(v));
a1 = x(2)*sin(w1*(1:length(v)) + x(4));
a2 = x(3)*sin(w2*(1:length(v)) + x(5));
p1 = x(4);
p2 = x(5);


t = (0:(length(v)-1))/Fs;
comps = nan([3, length(v)]);
comps(1, :) = d;
comps(2, :) = a1.*sin(w1*t + p1);
comps(3, :) = a2.*sin(w2*t + p2);

opt.r = 100;
opt.a = [a1(:)';a2(:)'];
opt.p = [p1*ones([1,length(v)]);p2*ones([1,length(v)])];
end