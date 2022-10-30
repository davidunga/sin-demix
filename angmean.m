function mu = angmean(angs)
% average angle (radians)
mu = angle(sum(exp(1i*angs)));