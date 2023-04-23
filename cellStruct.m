function s = cellStruct(params)
% Make struct of cell parameters

arguments
    params.c    % capacitance
    params.ve   % excitatory reversal potential
    params.vi   % inhibitory reversal potential
end

s = struct(c=params.c, ve=params.ve, vi=params.vi);
