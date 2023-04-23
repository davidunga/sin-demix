function res = stable_estimate_conductances(I,V,rest_win,expr,cell)


for i = 1 : 1
    disp(i);
    f1 = expr.f1 + .5 * 2 * (rand(1) - .5);
    f2 = expr.f2 + .5 * 2 * (rand(1) - .5);
    %f1 = expr.f1 + (pmin + rand(1) * (pmax - pmin));
    %f2 = expr.f2 + (pmin + rand(1) * (pmax - pmin));
    r = estimate_conductances(I,V,rest_win,exprStruct(Fs=expr.Fs,f1=f1,f2=f2),cell);
    if i == 1
        s = r;
        res = s;
    else
        flds = fieldnames(r);
        for k = 1 : length(flds)
            fld = flds{k};
            if isstruct(r.(fld))
                continue;
            end
            s.(fld) = s.(fld) + r.(fld);
            res.(fld) = s.(fld) / i;
        end
    end
end

