function comps = mixparams2comps(mp)

comps = nan([length(mp.ws)+1, length(mp.t)]);
comps(1,:) = mp.dc;
for i = 1 : length(mp.ws)
    comps(i+1, :) = mp.a(i,:).*sin(mp.ws(i)*mp.t + mp.p(i,:));
end
assert(~any(isnan(comps(:))));