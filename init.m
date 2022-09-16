function init()
% Initialize enviroment

BLACKLIST = {'deprecate', [ filesep '_' ], [ filesep '.' ]};

fprintf('Initializing %s ... ', cd);
restoredefaultpath();
subdirs = textscan(genpath(cd),'%s','Delimiter',pathsep);
subdirs = subdirs{1};
inds = true;
for k = 1 : length(BLACKLIST)
    inds = inds & ~contains(subdirs,lower(BLACKLIST{k}));
end
addpath(sprintf('%s;',subdirs{inds}));

fprintf('Done.\n');