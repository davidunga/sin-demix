function [a,ii] = getmid(a,d)

ix = floor((length(a) - d)/2)+1;
ii = ix:ix+d-1;
a = a(ii);