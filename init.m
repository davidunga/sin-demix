function init()
% Initialize enviroment

IGNORE = ["deprecate", filesep + "_", filesep + "."];

fprintf('--- %s --- \n', upper(cd));
fprintf('Initializing paths... ');
restoredefaultpath();
subdirs = split(string(genpath(cd)),pathsep);
addpath(join(subdirs(~contains(subdirs,IGNORE)),pathsep));
fprintf('Done.\n');