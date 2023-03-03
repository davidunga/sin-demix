function s = cellStruct(params)
% make struct of cell parameters

arguments
    params.c double
    params.ve double
    params.vi double
end

s = struct(c=params.c, ve=params.ve, vi=params.vi);
