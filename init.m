function init()
% Initialize enviroment

set(0, 'DefaultFigureRenderer', 'painters');

IGNORE = ["deprecate", filesep + "_", filesep + "."];

fprintf('--- %s --- \n', upper(cd));
fprintf('Initializing paths... ');
restoredefaultpath();
subdirs = split(string(genpath(cd)),pathsep);
addpath(join(subdirs(~contains(subdirs,IGNORE)),pathsep));
fprintf('Done.\n');