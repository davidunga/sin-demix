function comps = params2comps(params, t, options)
% build components from params
% OUTPUT:
%   comps - [3, length(t)] components matrix: [dc; sin1; sin2]

arguments
    params                  % params struct
    t                       % time vec
    options.smooth = false  % apply smoothing?
    options.r = 0           % boundary margin size: replicate first and last [r] elements
end

t = t(:)';
for fld = string(fieldnames(params)')
    params.(fld) = expand_param(params.(fld), length(t), options.r);
end

if options.smooth
    Fs = (length(t)-1)/(t(end) - t(1));
    win = 1 * (2*pi)/max(params.w1(1),params.w2(1)) * Fs;
    params.a0 = smoothdata(params.a0,"movmean",win);
end

comps = [
    params.a0;
    params.a1 .* sin(params.w1 .* t + params.p1);
    params.a2 .* sin(params.w2 .* t + params.p2)
    ];

if options.r > 0
    r = options.r;
    comps(:,1:r) = repmat(comps(:,r+1),[1,r]);
    comps(:,end-(r-1):end) = repmat(comps(:,end-r),[1,r]);
end


function x = expand_param(x, L, r)
if length(x) == 1
    x = x * ones(1,L);
elseif r > 0
    x(1:r) = x(r+1);
    x(end-r:end) = x(end-r);
end
assert(length(x) == L);
