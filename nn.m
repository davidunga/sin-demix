

x = [ret.I; ret.V];
y = [ret.gt.GTOT; ret.gt.GE];

X = {};
Y = {};
for i = 1 : 1000
    X{end+1} = x + randn(size(x));
    Y{end+1} = y + randn(size(y));
end

layers = [
    sequenceInputLayer
    convolution1dLayer(100, 8, Padding='same')
    batchNormalizationLayer
    reluLayer
    convolution1dLayer(100, 2, Padding='same')
    batchNormalizationLayer
    reluLayer
    ];