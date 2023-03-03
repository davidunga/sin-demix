function s = exprStruct(params)
% make struct of experiment parameters

arguments
    params.Fs double
    params.f1 double = NaN
    params.f2 double = NaN   
end

assert(isnan(params.f1) == isnan(params.f2));

s = struct();
s.Fs = params.Fs;
s.w1 = params.f1 * 2 * pi;
s.w2 = params.f2 * 2 * pi;
